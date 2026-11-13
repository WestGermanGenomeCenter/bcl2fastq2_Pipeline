import os
import re
import sys
import pandas as pd
from snakemake.exceptions import WorkflowError

p = os.path.abspath(".")
sys.path.append(os.path.abspath(os.getcwd())+"/scripts")
from helper import validateSamplesheet, validateOutput, validateProjectNum

# config.yaml may already have been loaded by the including snakefile, loading it again is harmless
configfile: "config.yaml"
outputfolder = config["demux"]["OutputFolder"]
samplesheet = config["demux"]["SampleSheet"]   # needed below, defined here so this file does not depend on the including one

# Check for projectNum in samplesheet name
projectNum=validateProjectNum(samplesheet)
if not projectNum:
    print("There was no project number to be found in the samplesheet and/or the name of the samplesheet")

# ---------------------------------------------------------------------------------------
# Read layout (SE / PE) detection
#
# A read tag is only looked for at the END of the name: _R1 / _R2 (read) or _I1 / _I2 (index read),
# optionally followed by Illumina's _001. Examples: S1_L001_R1_001, S1_R2, S1_I1_001.
# An "_R" anywhere else in a sample name (_Rep1, _RNA, ...) is ignored.
# Not supported: bare _1/_2 suffixes (ambiguous with sample names like "Patient_1").
# ---------------------------------------------------------------------------------------
READ_TAG = re.compile(r"^(?P<base>.+)_(?P<kind>[RI])(?P<num>\d+)(?P<suffix>_001)?$")

def splitReadTag(name):
    """'X_S1_R1_001' -> ('X_S1', 'R', 1, '_001'); name without read tag -> (name, None, None, '')"""
    m = READ_TAG.match(name)
    if not m:
        return name, None, None, ""
    return m.group("base"), m.group("kind"), int(m.group("num")), m.group("suffix") or ""

# fastq names without ".fastq.gz"; index reads (_I1/_I2 tag) are dropped here, once, for the whole pipeline
sample_names = list()
index_read_names = list()
if os.path.isfile(outputfolder+"/fastq_infiles_list.tx"):
    samples_dataframe = pd.read_csv(outputfolder+"/fastq_infiles_list.tx", header=None)
    for fastq in samples_dataframe.iloc[:, 0].values:
        name = fastq[:-9] if fastq.endswith(".fastq.gz") else fastq
        if splitReadTag(name)[1] == "I":
            index_read_names.append(name)
        else:
            sample_names.append(name)
    if index_read_names:
        print("Detected", len(index_read_names), "index read file(s) (_I1/_I2), ignoring them.")


def _detect_single_end() -> bool:
    """
    True = single-end, False = paired-end. Evaluated once at parse time.
    Paired-end only if every fastq has an _R1/_R2 tag and every R1 has its R2 (and vice versa).
    If one sample is single-end or one R2 is missing, the whole set is treated as single-end (a warning lists the culprits).
    Note: sample_names comes from fastq_infiles_list.tx, which exists only after the demux pipeline ran.
    On the very first invocation the list is empty and this returns True (single-end).
    """
    reads = {1: dict(), 2: dict()}
    untagged, unexpected = list(), list()
    for name in sample_names:
        base, kind, num, suffix = splitReadTag(name)
        if kind is None:
            untagged.append(name)
        elif num in reads:
            reads[num][base] = name
        else:
            unexpected.append(name)

    if unexpected:
        raise WorkflowError("Unexpected read tag(s) (only _R1/_R2 are read, _I1/_I2 are ignored): "
                            + ", ".join(unexpected[:5]) + (" ..." if len(unexpected) > 5 else "")
                            + ". Rename/remove these files (e.g. a third read R3) or use skip_demux with R1/R2 only.")
    if not reads[1] and not reads[2]:
        return True                                   # nothing tagged (or empty list): single-end
    if untagged:
        print("WARNING: fastq(s) without _R1/_R2 tag next to tagged ones:", ", ".join(untagged[:5]),
              "-> treating the whole set as single-end.")
        return True
    unpaired = sorted(set(reads[1]) ^ set(reads[2]))
    if unpaired:
        print("WARNING: no matching R1/R2 mate for:", ", ".join(unpaired[:5]) + (" ..." if len(unpaired) > 5 else ""),
              "-> treating the whole set as single-end.")
        return True
    # paired-end: the PE rules address files as <base>_R1_001.fastq.gz / <base>_R2_001.fastq.gz
    bad = [n for r in (1, 2) for n in reads[r].values() if not n.endswith("_R%d_001" % r)]
    if bad:
        raise WorkflowError("Paired-end input detected, but the PE rules require Illumina naming (<sample>_R1_001.fastq.gz / "
                            "<sample>_R2_001.fastq.gz). Rename e.g.: " + ", ".join(bad[:3]) + (" ..." if len(bad) > 3 else ""))
    return False

_SINGLE_END = _detect_single_end()


def isSingleEnd() -> bool:
    """True = single-end, False = paired-end (decided once, see _detect_single_end)"""
    return _SINGLE_END


def getSample_names_post_mapping():
    """SE: the fastq names. PE: one name per pair, X_R1_001 -> X_pe_001 (matches the align_pe output names)."""
    if isSingleEnd():
        return sample_names
    return [base + "_pe_001" for base, kind, num, suffix in map(splitReadTag, sample_names) if num == 1]


def getFiles():
    files = list()
    for sample in sample_names:
        trySplit = sample.split(os.sep)[-1]
        if trySplit:
            files.append(splitReadTag(trySplit)[0])
    return files


def _fastq_paths(files):
    """Input fastq paths for the given sample names, depending on which preprocessing steps are active."""
    if config["umi_tools"]["umi_tools_active"]:
        return expand("{out}/umi_extract/{file}.umis-extracted.fastq.gz", out=outputfolder, file=files)
    if config["cutadapt"]["cutadapt_active"]:
        return expand("{out}/trimmed/{file}_trimmed.fastq.gz", out=outputfolder, file=files)
    return expand("{out}/untrimmed_fastq/{file}.fastq.gz", out=outputfolder, file=files)


def getSamples(wildcards):
    if isSingleEnd():
        w = str(wildcards)
        # exact match first, so sample "S1_R1_001" is not shadowed by "S1_R1_001x"; fall back to substring as before
        hits = [s for s in sample_names if s == w] or [s for s in sample_names if w in s]
        if not hits:
            raise WorkflowError("No fastq found for wildcard '%s'. Known samples: %s"
                                % (w, ", ".join(sample_names[:5]) + (" ..." if len(sample_names) > 5 else "")))
        return _fastq_paths(hits[0])
    # paired-end: all R1 and R2 files (the wildcard is not used to select a pair)
    return _fastq_paths(list(sample_names))


def getFastqs(wildcards):
    samples = getSamples(wildcards)
    if len(samples) > 1:
        return " ".join(samples)
    return samples