from enum import Enum

class VariableKind(Enum):
    Atom="Atom"
    Bond = "Bond"
    Nonbond="Nonbond"
    Clash="Clash"
    Angle = "Angle"

class VariableID:
    @staticmethod
    def Atom(chunk:'AtomChunk'):
        return VariableID(
            chunk.get_disordered_tag(), #f"{chunk.resnum}.{chunk.name}",
            VariableKind.Atom,
            chunk.get_site_num(),
            chunk.is_water,
            chunk.name)
    def __init__(self,name:str,kind:VariableKind,site_num=None,is_water=None,atom_name=None):
        self.name=str(name)
        self.atom_name=atom_name
        self.kind = kind
        self.site_num=site_num
        self.is_water = is_water
    def __repr__(self):
        return self.name
    def __hash__(self):
        return hash((self.name, self.kind))
    def __eq__(self,other):
        return (self.name,self.kind) == (other.name,other.kind) 
    def __ne__(self,other):
        return not(self == other) 