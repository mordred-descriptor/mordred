"""Acid Base descriptor.

References:
    * http://cdk.github.io/cdk/1.5/docs/api/org/openscience/cdk/qsar/descriptors/molecular/AcidicGroupCountDescriptor.html
    * http://cdk.github.io/cdk/1.5/docs/api/org/openscience/cdk/qsar/descriptors/molecular/BasicGroupCountDescriptor.html

"""  # noqa: E501

from abc import abstractproperty

from rdkit import Chem

from ._base import Descriptor

__all__ = ("AcidicGroupCount", "BasicGroupCount")


class SmartsCountBase(Descriptor):
    __slots__ = ("_mol",)

    @classmethod
    def preset(cls, version):
        yield cls()

    def _create_smarts(self):
        s = ",".join("$(" + s + ")" for s in self.SMARTS)
        self._mol = Chem.MolFromSmarts("[" + s + "]")
        return self._mol

    @abstractproperty
    def SMARTS(self):
        pass

    def __str__(self):
        return self._name

    def parameters(self):
        return ()

    def calculate(self):
        pat = getattr(self, "_mol", None) or self._create_smarts()
        return len(self.mol.GetSubstructMatches(pat))

    rtype = int

    _extra_docs = ("SMARTS",)


class AcidicGroupCount(SmartsCountBase):
    r"""acidic group count descriptor."""

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "acidic group count"

    _name = "nAcid"

    SMARTS = (
        "[O;H1]-[C,S,P]=O",
        "[*;-;!$(*~[*;+])]",
        "[NH](S(=O)=O)C(F)(F)F",
        "n1nnnc1",
    )


class BasicGroupCount(SmartsCountBase):
    r"""basic group count descriptor."""

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "basic group count"

    _name = "nBase"

    SMARTS = (
        "[NH2]-[CX4]",
        "[NH](-[CX4])-[CX4]",
        "N(-[CX4])(-[CX4])-[CX4]",
        "[*;+;!$(*~[*;-])]",
        "N=C-N",
        "N-C=N",
    )
