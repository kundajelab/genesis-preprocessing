// Checks modules/common/retry.nf offline for tests/verify_nextflow_style.py. With --mode functions,
// it reports each function's result on the task logs under --logs (one folder per case); with
// --mode oom or --mode other, a task fails its first attempt with exit 1, after Slurm's oom_kill
// message (oom) or without it (other), and succeeds on any retry.
include { slurmReportedOutOfMemory ; wasOutOfMemory ; oomMemoryOf ; oomErrorStrategy } from '../modules/common/retry'

process ATTEMPT {
    memory { oomMemoryOf(task, 100.MB) }
    errorStrategy { oomErrorStrategy(task) }
    maxRetries 2
    input:
    val kind
    output:
    stdout
    script:
    """
    echo "attempt ${task.attempt} memory ${task.memory.toMega()}"
    if [ ${task.attempt} = 1 ]; then
        if [ '${kind}' = oom ]; then
            echo '[2026-11-08T04:24:28.700] error: Detected 1 oom_kill event in StepId=1.batch. Some of the step tasks have been OOM Killed.' >&2
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
     'tail-oom': true].each { name, want ->
        results << report("slurmReportedOutOfMemory(${name})", slurmReportedOutOfMemory("${logs}/${name}"), want)
    }
    results << report('slurmReportedOutOfMemory(null)', slurmReportedOutOfMemory(null), false)
    [[1, 'oom', true], [1, 'other', false], [137, 'other', true], [135, 'other', true],
     [140, 'oom', true], [140, 'other', false], [143, 'oom', false], [143, 'other', false],
     [0, 'oom', false]].each { status, dir, want ->
        results << report("wasOutOfMemory(${status}, ${dir})", wasOutOfMemory(status, "${logs}/${dir}"), want)
    }
    def afterOom = [previousTrace: [memory: 8.GB.toBytes(), workdir: "${logs}/oom"], exitStatus: 1]
    def afterOther = [previousTrace: [memory: 8.GB.toBytes(), workdir: "${logs}/other"], exitStatus: 143]
    results << report('memory, first attempt', oomMemoryOf([previousTrace: null], 17.GB), 17.GB)
    results << report('memory after running out', oomMemoryOf(afterOom, 17.GB), 16.GB)
    results << report('memory after another failure', oomMemoryOf(afterOther, 17.GB), 8.GB)
    def atCeiling = [exitStatus: 1, workDir: "${logs}/oom", memory: 247.GB, resourceLimits: [memory: 247.GB]]
    results << report('retry at the memory ceiling', oomErrorStrategy(atCeiling), 'terminate')
    println "RESULT ${results.count { r -> !r }} wrong of ${results.size()}"
}
