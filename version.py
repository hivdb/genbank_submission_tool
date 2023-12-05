import sys
from packaging import version

MIN_VERSION = "3.9"


VERSION = '.'.join(map(str, sys.version_info[:3]))

if version.parse(VERSION) >= version.parse(MIN_VERSION):
    sys.exit(0)
else:
    print(f'Version is {VERSION}, lower than {MIN_VERSION}.')
    sys.exit(1)
