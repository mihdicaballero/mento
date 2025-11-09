from importlib.metadata import version, PackageNotFoundError

def get_version(dist_name: str = "mento") -> str:
    try:
        return version(dist_name)
    except PackageNotFoundError:
        return "0+local"

__version__ = get_version()