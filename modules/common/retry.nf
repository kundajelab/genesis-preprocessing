// Out-of-memory retries for processes whose memory needs are hard to predict. Such a process
// estimates its memory from its inputs and retries only when an attempt ran out of memory, with
// twice the memory each time; any other failure stops the run as before:
//
//     include { oomMemoryOf ; oomErrorStrategy } from '../common/retry'
//
//     process EXAMPLE {
//         memory { oomMemoryOf(task, 2.GB + 4 * bam.size().B) }
//         errorStrategy { oomErrorStrategy(task) }
//         maxRetries 2
//
// NOTE: these must be called inside a directive closure. A bare call in a process body is taken
// for a directive name and fails with "Unknown process directive".

// Whether Slurm recorded killing processes of the job in workDir for running out of memory, in
// the job's .command.log: "error: Detected N oom_kill event(s) in StepId=...". Only the end of
// the log is read, which is where Slurm writes it.
def slurmReportedOutOfMemory(workDir) {
    if (workDir == null) {
        return false
    }
    def logFile = java.nio.file.Paths.get(workDir.toString(), '.command.log')
    try {
        if (!java.nio.file.Files.exists(logFile)) {
            return false
        }
        // NOTE: a ByteBuffer rather than a byte array, which Nextflow's script syntax cannot create
        java.nio.channels.FileChannel.open(logFile).withCloseable { log ->
            def tailSize = Math.min(log.size(), 65536L)
            def tail = java.nio.ByteBuffer.allocate(tailSize as int)
            // A read of a regular file fills the buffer, up to the end of the file.
            log.read(tail, log.size() - tailSize)
            return new String(tail.array(), 0, tail.position(), 'UTF-8').contains('oom_kill')
        }
    } catch (IOException _exception) {
        return false
    }
}

// Whether an attempt ran out of memory: killed with SIGKILL (137), or SIGBUS (135), which a job
// can get when it touches memory the kernel cannot provide. Or, when Slurm logged that it killed
// some of the job's processes for memory:
// - exit status 1: a script that handles a failed command (as QC does SPP's) or a killed trace
//   monitor then exits 1, hiding the 137;
// - exit status 140, the time limit's SIGUSR2 on Slurm: a program whose worker processes were
//   killed can wait for them until then, as bamCoverage does. Without an OOM kill, a 140 is
//   only out of time, which more memory would not help.
def wasOutOfMemory(exitStatus, workDir) {
    exitStatus in [135, 137] || (exitStatus in [1, 140] && slurmReportedOutOfMemory(workDir))
}

// Memory for this attempt: the estimate on the first, then twice the last attempt's when it ran
// out of memory, or the same when it failed otherwise.
def oomMemoryOf(task, baseMem) {
    if (!task.previousTrace) {
        return baseMem
    }
    def previous = task.previousTrace.memory as long
    def previousOom = wasOutOfMemory(task.exitStatus, task.previousTrace.get('workdir'))
    return (previousOom ? 2 * previous : previous).B
}

// Retry an attempt only when it ran out of memory, and not when it already had the most memory
// the profile allows (resourceLimits): the next attempt would get the same again.
def oomErrorStrategy(task) {
    def ceiling = task.resourceLimits?.memory
    def atCeiling = ceiling && task.memory && task.memory >= ceiling
    wasOutOfMemory(task.exitStatus, task.workDir) && !atCeiling ? 'retry' : 'terminate'
}
