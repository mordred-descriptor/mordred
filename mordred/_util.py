def parse_enum(enum, v):
    if isinstance(v, enum):
        return v
    else:
        return enum[v]
