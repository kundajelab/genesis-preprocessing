// Checks modules/common/retry.nf offline for tests/verify_nextflow_style.py. With --mode functions,
// it reports each function's result on the task logs under --logs (one folder per case); with
// --mode oom or --mode other, a task fails its first attempt with exit 1, after Slurm's oom_kill
// message (oom) or without it (other), and with --mode timeout with exit 140 after Slurm's time
// limit message; it succeeds on any retry.
include { slurmReportedOutOfMemory ; slurmReportedTimeLimit ; wasOutOfMemory ; wasOutOfTime ; oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../modules/common/retry'

process ATTEMPT {
    memory { oomMemoryOf(task, 100.MB) }
    time { retryTimeOf(task, 10.min) }
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    input:
    val kind
    output:
    stdout
    script:
    """
    echo "attempt ${task.attempt} memory ${task.memory.toMega()} time ${task.time.toMinutes()}"
    if [ ${task.attempt} = 1 ]; then
        if [ '${kind}' = oom ]; then
            echo '[2026-11-08T04:24:28.700] error: Detected 1 oom_kill event in StepId=1.batch. Some of the step tasks have been OOM Killed.' >&2
        elif [ '${kind}' = timeout ]; then
            echo 'slurmstepd: error: *** JOB 61234089 ON nid004212 CANCELLED AT 2026-11-13T12:44:25 DUE TO TIME LIMIT ***' >&2
            exit 140
        fi
        exit 1
    fi
    """
}

def report(label, got, want) {
    println "${got == want ? 'ok' : 'WRONG'} ${label}: ${got} (want ${want})"
    return got == want
}

workflow {
    if (params.mode != 'functions') {
        ATTEMPT(channel.of(params.mode)).view()
    } else {
        checkFunctions(params.logs)
    }
}

def checkFunctions(logs) {
    def results = []
    [oom: true, other: false, missing: false, empty: false, 'head-only': false,
     'tail-oom': true, timeout: false].each { name, want ->
        results << report("slurmReportedOutOfMemory(${name})", slurmReportedOutOfMemory("${logs}/${name}"), want)
    }
    results << report('slurmReportedOutOfMemory(null)', slurmReportedOutOfMemory(null), false)
    [timeout: true, oom: false, other: false, missing: false].each { name, want ->
        results << report("slurmReportedTimeLimit(${name})", slurmReportedTimeLimit("${logs}/${name}"), want)
    }
    results << report('slurmReportedTimeLimit(null)', slurmReportedTimeLimit(null), false)
    [[1, 'oom', true], [1, 'other', false], [137, 'other', true], [135, 'other', true],
     [140, 'oom', true], [140, 'other', false], [143, 'oom', false], [143, 'other', false],
     [0, 'oom', false]].each { status, dir, want ->
        results << report("wasOutOfMemory(${status}, ${dir})", wasOutOfMemory(status, "${logs}/${dir}"), want)
    }
    // A job killed before it recorded its exit status has none (Integer.MAX_VALUE in Nextflow).
    def unknown = Integer.MAX_VALUE
    [[140, 'timeout', true], [140, 'other', true], [140, 'oom', false], [unknown, 'timeout', true],
     [unknown, 'other', false], [1, 'other', false], [137, 'timeout', false],
     [0, 'other', false]].each { status, dir, want ->
        results << report("wasOutOfTime(${status}, ${dir})", wasOutOfTime(status, "${logs}/${dir}"), want)
    }
    def afterOom = [previousTrace: [memory: 8.GB.toBytes(), workdir: "${logs}/oom"], exitStatus: 1]
    def afterOther = [previousTrace: [memory: 8.GB.toBytes(), workdir: "${logs}/other"], exitStatus: 143]
    results << report('memory, first attempt', oomMemoryOf([previousTrace: null], 17.GB), 17.GB)
    results << report('memory after running out', oomMemoryOf(afterOom, 17.GB), 16.GB)
    results << report('memory after another failure', oomMemoryOf(afterOther, 17.GB), 8.GB)
    def afterTimeout = [previousTrace: [time: 3.h.toMillis(), memory: 8.GB.toBytes(), workdir: "${logs}/timeout"], exitStatus: 140]
    def afterOomTime = [previousTrace: [time: 3.h.toMillis(), workdir: "${logs}/oom"], exitStatus: 1]
    results << report('time, first attempt', retryTimeOf([previousTrace: null], 5.h), 5.h)
    results << report('time after running out', retryTimeOf(afterTimeout, 5.h), 6.h)
    results << report('time after running out of memory', retryTimeOf(afterOomTime, 5.h), 3.h)
    results << report('memory after running out of time', oomMemoryOf(afterTimeout, 17.GB), 8.GB)
    def atCeiling = [exitStatus: 1, workDir: "${logs}/oom", memory: 247.GB, resourceLimits: [memory: 247.GB]]
    results << report('retry at the memory ceiling', retryErrorStrategy(atCeiling), 'terminate')
    def timedOut = [exitStatus: 140, workDir: "${logs}/timeout", memory: 247.GB, time: 24.h,
                    resourceLimits: [memory: 247.GB, time: 48.h]]
    results << report('retry out of time', retryErrorStrategy(timedOut), 'retry')
    def timeAtCeiling = timedOut + [time: 48.h]
    results << report('retry at the time ceiling', retryErrorStrategy(timeAtCeiling), 'terminate')
    def failed = [exitStatus: 1, workDir: "${logs}/other", memory: 1.GB, time: 1.h, resourceLimits: [:]]
    results << report('retry another failure', retryErrorStrategy(failed), 'terminate')
    println "RESULT ${results.count { r -> !r }} wrong of ${results.size()}"
}
