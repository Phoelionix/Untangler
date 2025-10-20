from untangle import Untangler


untangler = Untangler()

#untangler.prepare_pdb_and_read_altlocs("output/refmacout_minRfree_loopEnd2.pdb","refmacout_minRfree_loopEnd2_sepchain.pdb",sep_chain_format=True)
untangler.prepare_pdb_and_read_altlocs("data/HoltonbestRfree.pdb","HoltonbestRfreeOut.pdb",sep_chain_format=True,altloc_from_chain_fix=True)