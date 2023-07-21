from packaging.version import Version as StrictVersion

from mordred.RingCount import RingCount

# Start Code 4
presets = list(RingCount.preset(version=StrictVersion("1.0.0")))
print(len(presets))
