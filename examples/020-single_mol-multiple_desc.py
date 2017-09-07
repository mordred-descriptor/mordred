from rdkit import Chem

from mordred import Chi, ABCIndex, RingCount, Calculator, is_missing, descriptors

benzene = Chem.MolFromSmiles("c1ccccc1")

# Create empty Calculator instance
calc1 = Calculator()

# Register descriptor instance
calc1.register(Chi.Chi(type="path_cluster", order=4))

# Register descriptor class using preset
calc1.register(RingCount.RingCount)

# Register all descriptors in module
calc1.register(ABCIndex)


# Calculate descriptors
result = calc1(benzene)

print(result)
# >>> [0.0, 1, 0, 0, 0, 1, (snip)


# Calculator constructor can register descriptors
calc2 = Calculator(Chi.Chi)

# Descriptors module contains all descriptors
calc3 = Calculator(descriptors)

# User can access all descriptor instances by descriptors property
print(calc3.descriptors)
# >>> (mordred.EccentricConnectivityIndex.EccentricConnectivityIndex(), (snip)


# Calculate descriptors
result = calc3(benzene)

# get first missing value
na1 = next(r for r in result if is_missing(r))

# get reason
print(na1.error)
# >>> missing 3D coordinate


# Delete all missing value
result = result.drop_missing()


# convert to dict
print(result.asdict())
