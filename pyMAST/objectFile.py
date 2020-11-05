from anndata import AnnData;

class MAST:
    """
        Creates a MAST Object from scanpy object to work upon
    """
    def __init__(self, sc: AnnData):
        self.data = sc.X
        self.var_names = sc.var_names
        self.obj_names = sc.obs_names
    
    def __str__(self):
        return f"type('MAST') Object"
        
        