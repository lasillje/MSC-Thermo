def write_to_log(filename, message: str, mode="a", skip=True):
    if not skip:
        if filename is not None or filename != "" or filename:
            with open(filename, mode) as out:
                out.write(f"{message}\n")