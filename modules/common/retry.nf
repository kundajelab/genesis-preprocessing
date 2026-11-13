// Retries for processes whose memory and time needs are hard to predict. Such a process
// estimates both from its inputs and retries only when an attempt ran out of memory or of time,
// with twice that resource each time; any other failure stops the run as before:
//
//     include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'
//
//     process EXAMPLE {
//         memory { oomMemoryOf(task, 2.GB + 4 * bam.size().B) }
//         time { retryTimeOf(task, 10.min + 1.h * bam.size() / 1e9) }
//         errorStrategy { retryErrorStrategy(task) }
//         maxRetries 2
//
// NOTE: these must be called inside a directive closure. A bare call in a process body is taken
// for a directive name and fails with "Unknown process directive".

// Whether the end of the Slurm log of the job in workDir, its .command.log, contains text. Only
// the end is read, which is where Slurm writes its reports.
def slurmLogEndContains(workDir, text) {
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
            return new String(tail.array(), 0, tail.position(), 'UTF-8').contains(text)
        }
    } catch (IOException _exception) {
        return false
    }
}

// Whether Slurm recorded killing processes of the job in workDir for running out of memory:
// "error: Detected N oom_kill event(s) in StepId=...".
def slurmReportedOutOfMemory(workDir) {
    slurmLogEndContains(workDir, 'oom_kill')
}

// Whether Slurm recorded cancelling the job in workDir at its time limit:
// "*** JOB N ON NODE CANCELLED AT TIME DUE TO TIME LIMIT ***".
def slurmReportedTimeLimit(workDir) {
    slurmLogEndContains(workDir, 'DUE TO TIME LIMIT')
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

// Whether an attempt ran out of time. Nextflow has Slurm send the job SIGUSR2 30 s before its
// time limit (--signal B:USR2@30), on which the job script exits 140, as the DAP-seq run's timed
// out tasks did. A job still running at the limit is cancelled before it records an exit status,
// with Slurm's report in its log. Not when it also ran out of memory (see wasOutOfMemory): a
// program can wait for workers killed for memory until the time limit, and more time would not
// help it.
def wasOutOfTime(exitStatus, workDir) {
    !wasOutOfMemory(exitStatus, workDir) && (exitStatus == 140 || slurmReportedTimeLimit(workDir))
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

// Time for this attempt: the estimate on the first, then twice the last attempt's when it ran
// out of time, or the same when it failed otherwise.
def retryTimeOf(task, baseTime) {
    if (!task.previousTrace) {
        return baseTime
    }
    def previous = task.previousTrace.time as long
    def previousTimeout = wasOutOfTime(task.exitStatus, task.previousTrace.get('workdir'))
    return (previousTimeout ? 2 * previous : previous).ms
}

// Retry an attempt only when it ran out of memory or of time, and not when it already had the
// most of that the profile allows (resourceLimits): the next attempt would get the same again.
def retryErrorStrategy(task) {
    def memoryCeiling = task.resourceLimits?.memory
    def timeCeiling = task.resourceLimits?.time
    def memoryAtCeiling = memoryCeiling && task.memory && task.memory >= memoryCeiling
    def timeAtCeiling = timeCeiling && task.time && task.time >= timeCeiling
    def retryMemory = wasOutOfMemory(task.exitStatus, task.workDir) && !memoryAtCeiling
    def retryTime = wasOutOfTime(task.exitStatus, task.workDir) && !timeAtCeiling
    retryMemory || retryTime ? 'retry' : 'terminate'
}
