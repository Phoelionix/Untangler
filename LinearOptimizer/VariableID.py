from enum import Enum

class VariableKind(Enum):
    Atom="Atom"
    Bond = "Bond"
    Nonbond="Nonbond"
    Clash="Clash"
    Angle = "Angle"
    Penalty = "Penalty"

class VariableID:
    @staticmethod
    def Atom(chunk:'AtomChunk'):
        return VariableID(
            chunk.get_disordered_tag(), #f"{chunk.resnum}.{chunk.name}",
            VariableKind.Atom,
            chunk.get_site_num(),
            chunk.is_water,
            chunk.name,
            chunk.get_resname())
    def __init__(self,name:str,kind:VariableKind,site_num=None,is_water=None,atom_name=None,residue_name=None):
        self.name=str(name)
        self.kind = kind
        self.site_num=site_num
        self.is_water = is_water
        self.atom_name=atom_name
        self.resname=residue_name
    def get_site_num(self):
        return self.site_num
    def get_resname(self):
        return self.resname
    ###
    def __repr__(self):
        return self.name
    def __hash__(self):
        return hash((self.name, self.kind))
    def __eq__(self,other):
        return (self.name,self.kind) == (other.name,other.kind) 
    def __ne__(self,other):
        return not(self == other) 