def write_to_log(filename: str, message: str, mode="a"):
    with open(filename, mode) as out:
        out.write(f"{message}\n")