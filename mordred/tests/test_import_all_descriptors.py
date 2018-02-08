from mordred import Calculator, descriptors


def test_import_all_descriptors():
    calc = Calculator(descriptors)
    it = iter(calc.descriptors)
    before = next(it).__module__
    for current in it:
        current = current.__module__
        assert before <= current, "{!r} > {!r}".format(before, current)

        before = current
