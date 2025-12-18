import UntangleFunctions 
from Bio.PDB.Atom import Atom,DisorderedAtom

class DisorderedTag():
    def __init__(self,res_num,atom_name):
        self._resnum, self._name = int(res_num),str(atom_name)
    def is_entry(self,pdb_entry:UntangleFunctions.PDB_Atom_Entry)->bool:
        if not pdb_entry.valid:
            return False
        return DisorderedTag(pdb_entry.res_num,pdb_entry.atom_name)==self
    def is_riding_H_entry(self,H_pdb_entry:UntangleFunctions.PDB_Atom_Entry)->bool:
        if not H_pdb_entry.valid:
            return False
        if not H_pdb_entry.atom_name[0]=="H":
            return False
        if not int(H_pdb_entry.res_num)==self.resnum():
            return False
        name_length = len(self.atom_name())
        assert 1 <= name_length <= 4 
        if name_length<4:
            full_name = ' '+self.atom_name() + ' '*(3-name_length)
        else:
            full_name=self.atom_name()
        if UntangleFunctions.H_get_parent_fullname(H_pdb_entry.atom_name,[full_name],debug_none_return=False)==full_name:
            return True
        return False
    def resnum(self):
        return self._resnum
    def atom_name(self):
        return self._name 
    def element(self):
        return self.atom_name()[0]   # FIXME !!!!!!!!!!
    def num_bound_e(self):
        return UntangleFunctions.NUM_E[self.element()]
    def __repr__(self):
        return f"{self.resnum()}.{self.atom_name()}"
    # def __format__(self,format_spec):
    #     return str(self)
    # def __str__(self):
    #     return f"{self.resnum()}.{self.atom_name()}"
    
    def __hash__(self):
        return hash((self._resnum, self.atom_name()))

    def ordered_tag(self,altloc:str):
        return OrderedTag(self.resnum(),self.atom_name(),altloc)

    def __eq__(self, other:'DisorderedTag'):
        return (self._resnum, self.atom_name()) == (other._resnum, other.atom_name())
    def __ne__(self, other):
        return not(self == other)
    @staticmethod 
    def from_atom(a:Atom):
        def atom_res_seq_num(atom:Atom)->int:
                return atom.get_parent().get_id()[1]
        return DisorderedTag(atom_res_seq_num(a),a.get_name())

class OrderedTag(DisorderedTag):
    def __init__(self,res_num,atom_name,altloc):
        super().__init__(res_num,atom_name)
        self._altloc=altloc
    def altloc(self):
        return self._altloc
    def disordered_tag(self):
        return DisorderedTag(self.resnum(),self.atom_name())
    def __repr__(self):
        return f"{self.resnum()}.{self.atom_name()}.{self.altloc()}"
    def __hash__(self):
        return hash((self._resnum, self.atom_name(),self.altloc()))

    def __eq__(self, other:'OrderedTag'):
        return (self._resnum, self.atom_name(),self.altloc()) == (other._resnum, other.atom_name(),other.altloc())
    def __ne__(self, other):
        return not(self == other)
    @staticmethod 
    def from_atom(a:Atom):
        def atom_res_seq_num(atom:Atom)->int:
                return atom.get_parent().get_id()[1]
        return OrderedTag(atom_res_seq_num(a),a.get_name(),a.get_altloc())