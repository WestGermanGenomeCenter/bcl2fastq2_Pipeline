import os
import re


def validateSamplesheet(samplesheet):
    """
    Check for a few possible input errors and give an appropriate response
    """
    if not os.path.isfile(samplesheet):
        print("Couldn't be started. The following path or file %s does not exist" % (samplesheet,))
        return False
    return True


def validateOutput(outputfolder):
    """
    Check for a few possible input errors and give an appropriate response
    """
    if outputfolder.endswith('/'):
        outputfolder = outputfolder[:-1]
        if not os.access(dirname(outputfolder), os.W_OK):
            print("Couldn't be started. You don't have write permission for the following path %s" % (dirname(outputfolder),))
            print("It's also possible that parent folders/directories of the given path do not exist yet and need to be created")
            return None
    return outputfolder


def validateProjectNum(samplesheet):
    """
    Check if there is a projectNum in the samplesheet name
    """
    samplesheetName = os.path.basename(samplesheet)
    if "." in samplesheetName:
        samplesheetName = samplesheetName[:samplesheetName.index(".")]
    projectNum = re.search(r'\d+', samplesheetName)
    if projectNum:
        return projectNum.group(0)
    else:
        return ""


def adaptersToStringParams(adapters,type):
    """
    Convert list of adapters to string in form of concated params for cutadapt
    """
    if len(type)==len(adapters):
        adapter_params = ['-{} {}' for adapter in adapters]
        string = ' '.join('-{} {}'.format(*t) for t in zip(type,adapters))
        return string
    else:
        adapter_params = ['-'+type[0]+' {}' for adapter in adapters]
        string = ' '.join(adapter_params)
        return string.format(*adapters)

class SamplesHaveAWrongFormat(Exception):
    """Exception raised for errors in determining the output that should be moved."""
    pass
