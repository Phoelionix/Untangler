from LinearOptimizer.Input import ConstraintsHandler,DisorderedTag,OrderedAtomLookup
#from LinearOptimizer.Solver import VariableID,VariableKind
import UntangleFunctions
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom
import numpy as np


# Tension between Geo and X-ray


def create_delta_constraints_file(pre_unrestrained_constraints_file,post_unrestrained_constraints_file):
    pass


class GeoXrayTension:
    # Computes the deltas for geometries and atom positions. 
    # Tension is ratio of geometry deviation over sum of involved atom deltas # TODO Linear may not be best.  
    def __init__(self,pdb_files,symmetries,water_water_nonbond):
        print(f"Calculating tension reflected in changes between pdb files: {pdb_files}")
        assert len(set(pdb_files))==len(pdb_files)

        def separation(v1,v2):
            return np.sqrt(np.sum((v1-v2)**2))
        atom_residual_dicts:list[dict[DisorderedTag,float]]=[]
        atom_pos_dicts:list[dict[DisorderedTag,dict[str,tuple[float]]]]=[]
        for pdb in pdb_files:
            atom_pos_dict:dict[DisorderedTag,dict[str,tuple[float]]]={}
            constraints_handler = ConstraintsHandler()

            geo_file = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{UntangleFunctions.model_handle(pdb)}.geo"

            struct = PDBParser().get_structure("struct",pdb)
            ordered_atom_lookup=OrderedAtomLookup(struct.get_atoms(),waters=True)
            
            constraints_handler.load_all_constraints(geo_file,[geo_file],[],ordered_atom_lookup, symmetries,water_water_nonbond)
            atom_residual_dicts.append(constraints_handler.atom_residuals)


            for atom in ordered_atom_lookup.select_atoms_by():
                #disordered_tag = DisorderedTag(DisorderedTag(OrderedAtomLookup.atom_res_seq_num(atom),atom.get_name()),VariableKind.Atom)
                disordered_tag = DisorderedTag(OrderedAtomLookup.atom_res_seq_num(atom),atom.get_name())
                if disordered_tag not in atom_pos_dict:
                    atom_pos_dict[disordered_tag]={}
                atom_pos_dict[disordered_tag][atom.get_altloc()]=atom.get_coord()
            atom_pos_dicts.append(atom_pos_dict)

            
        # TODO tension should be calculated on a per conformer basis. Then average taken.
        print("Computing tensions")
        self.site_tensions:dict[DisorderedTag,float]={}
        assert len(atom_residual_dicts)>1
        assert len(atom_residual_dicts)==len(atom_pos_dicts)
        for key in atom_residual_dicts[0]:
            broke=False
            # skip sites (waters) that were removed # XXX better implementation?
            for d in atom_residual_dicts:
                if key not in d:
                    broke=True
                    break
            if broke:
                self.site_tensions[key]=0
                continue
            residuals = [d[key] for d in atom_residual_dicts]
            # for d in atom_pos_dicts:
            #     assert key in d, list(d.keys())
            positions = [d[key] for d in atom_pos_dicts]
            deltas_geo = [residuals[i+1]-residuals[i] for i in range(len(residuals)-1)]
            displacements = [] 
            for i in range(len(positions)-1):
                assert len(positions[i])==len(positions[i+1])
                mean_separation=0
                for altloc in positions[i]:
                    if altloc not in positions[i+1]:
                        continue  # TODO why is this necessary? Waters?
                    mean_separation += separation(positions[i+1][altloc],positions[i][altloc])
                mean_separation/=len(positions[i])
                displacements.append(mean_separation)
            self.site_tensions[key]=np.mean([GeoXrayTension.tension(g,d) for g,d in zip (deltas_geo, displacements)]) 
        mean_tension = np.mean([abs(t) for t in self.site_tensions.values()])
        assert mean_tension > 0
        for key in self.site_tensions:
            self.site_tensions[key]/=mean_tension
    @staticmethod
    # TODO delta_pos should be relative to average mean shift of all atoms involved in the geometry
    def tension(delta_geo,delta_pos):
        if delta_pos == 0:
            #return 1
            return 0
        #return max(0,delta_geo/delta_pos)
        return delta_geo/delta_pos

