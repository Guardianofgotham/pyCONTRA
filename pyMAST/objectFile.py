from anndata import AnnData;

class MAST:
    def __init__(self, sc: AnnData):
        self.data = sc.X
        
        