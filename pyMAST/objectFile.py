from anndata import AnnData;

class MAST:
    """
        Creates a MAST Object from scanpy object to work upon
    """
    def __init__(self, sc: AnnData):
        self.data = sc.X
        
        