import subprocess

def run_cmd(cmd):
    """    
    Arguments:
        cmd {list} -- Command to run.
    """
    logger.info("Running: %", " ".join(cmd))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
