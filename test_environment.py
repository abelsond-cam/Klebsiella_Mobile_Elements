"""Check that the current Python is 3.x. Used by 'make requirements' / test_environment (not a config file)."""
import sys


def main():
    system_major = sys.version_info.major
    required_major = 3

    if system_major != required_major:
        raise TypeError(
            "This project requires Python {}. Found: Python {}".format(
                required_major, sys.version
            )
        )
    else:
        print(">>> Development environment passes all tests!")


if __name__ == "__main__":
    main()
