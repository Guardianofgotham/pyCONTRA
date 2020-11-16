import os.path


def checkFileExists(filename: str):
    if(not os.path.exists(filename)):
        error(f"FileNotExists: {filename} not found")

def error(msg: str):
    raise Exception(msg)