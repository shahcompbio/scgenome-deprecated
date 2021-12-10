import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import anndata as ad

import seaborn
import numpy as np
import pandas as pd
import pylab
import sklearn.preprocessing
import scipy

import pyranges
import copy

import scgenome
import scgenome.utils
import scgenome.cnfilter


class NotDLPAnnDataType(Exception):
    def __init__(self, reason, message="Not DLP AnnData Type!"):
        self.reason = reason
        self.message = message
        super().__init__(f"{self.message} {self.reason}")

class DLPAnnData(ad.AnnData):
    def __init__(self):
        pass
    
    def __getitem__(self, index) -> "DLPAnnData":
        ndad = DLPAnnData()
        ndad.from_anndata(super().__getitem__(index))
        return ndad
    
    def write_h5ad(self, path, *args, **kwargs):
        """
        Writes current DLPAnnData into the provided path as AnnData h5ad.
        
        Parameters
        ------
        path
            Path to write current DLPAnnData into.
        """
        oldVarIndex = self.var.index.names
        oldObsIndex = self.obs.index.names
        
        self.var.reset_index(inplace=True)
        self.obs.reset_index(inplace=True)
        
        super().write(path, *args, **kwargs)
        
        self.var.set_index(oldVarIndex, inplace=True)
        self.obs.set_index(oldObsIndex, inplace=True)
        
    
    def load_h5ad(self, path, *args, **kwargs):
        """
        Loads h5ad from path into the current DLPAnnData.
        
        Parameters
        ------
        path
            Path to the h5ad file to load into current DLPAnnData.
        
        Returns
        ------
        None
            Instantiates DLPAnnData with default parameters infered from AnnData.
        """
        super().read_h5ad(path, *args, **kwargs)
        self.var.set_index(['chr', 'start', 'end'])
        
    def generate_pyranges(self):
        """
        Generates pyranges object in .uns field
        """
        df = self.to_df()

        self.uns = {}
        self.uns['pyranges'] = pyranges.from_dict({
                "Chromosome": df.columns.get_level_values('chr'),
                "Start": df.columns.get_level_values('start'),
                "End": df.columns.get_level_values('end'),
                "cell_id": self.obs.index.to_list(),
                "copy": self.layers['copy'].flatten(),
                "gc": self.layers['gc'].flatten(),
                "reads": self.layers['reads'].flatten(),
            })
        
    def set_X(self, property_to_set = 'counts', inplace = False):
        """
        Make X to the property that we have stored in `layers` and track that it is in X.

        Parameters
        ----------
        property_to_set
            The property to track inside the anndata object. Places in X and marks property is followed.
        inplace
            Return the new DLP_AnnData or write in place.
            
        Returns
        ------
        New DLPAnnData with X set to `property_to_set` as the X if not inplace.
        Otherwise, modifies self with X.
        
        Properties Set
        ------
        DLPAnnData.X - X array from layers
        DLPAnndata.current_X - the `property_to_set`
        """
        if property_to_set in self.layers.keys():
            if inplace:
                self.X = self.layers[property_to_set]
                self.current_X = property_to_set
            else:
                dadcopy = self.copy()
                dadcopy.set_X(property_to_set, inplace=True)
                return dadcopy
        else:
            print(f"WARNING: {property_to_set} is not in layers! Doing nothing.")
        
    def copy(self):
        """
        Makes a deep copy of the current DLPAnnData Object and returns.
        
        Returns
        ------
        Return deep copy DLP AnnData Object.
        """
        return copy.deepcopy( self )
        
        
    def from_anndata(self, adata):
        """
        Load DLP AnnData from an AnnData Object.
        Requires `reads`, `copy`, `state`.
        Optional `gc`.
        
        Parameters
        ------
        adata
            AnnData Object to wrap with DLP.
        
        Returns
        ------
        None
            Instantiates current DLPAnnData from the adata object.
        """
        
        required_layers = ['reads', 'copy', 'state']
        
        # Enforce required layers
        for layer in required_layers:
            if layer not in adata.layers.keys():
                raise NotDLPAnnDataType(f"{layer} not found in layers!")
            
        super().__init__(
            X= adata.X,
            obs=adata.obs,
            var=adata.var,
            layers=adata.layers,
        )
        
        # Check on optional layers
        if 'gc' not in adata.layers.keys():
            self.layers['gc'] = np.full(self.X.shape, 0.5)
                
        self.set_X('copy', inplace=True)
        
    def load_dlp_files(self, cn_file, metrics_file):
        """
        Load DLP from given directory paths.
        
        Parameters
        ------
        cn_dir
            cn data file
        metrics_dir
            
            
        Returns
        ------
        DLPAnnData
            AnnData wrapper for dealing with DLP Data.
        """
        cn_data = pd.read_csv(
            cn_file,
            dtype={
                'cell_id': 'category',
                'sample_id': 'category',
                'library_id': 'category',
                'chr': 'category',
            })

        metrics_data = pd.read_csv(
            metrics_file,
            dtype={
                'cell_id': 'category',
                'sample_id': 'category',
                'library_id': 'category',
            })

        self.load_dlp(cn_data, metrics_data)
    
    def load_dlp(self, cn, metrics):
        """
        Load DLPAnnData given the CN DataFrame and Metrics DataFrame.
        
        Parameters
        ------
        cn
            The CN DataFrame. This is required to have the fields `[chr, start, end, cell_id, reads, copy, state]`.
        metrics
            The metrics DataFrame. This is required to have `cell_id`. Becomes the obs of the AnnData.
            
        Returns
        ------
        None
            Instantiates current object
        """
        cn_data = cn
        metrics_data = metrics
        
        
        if 'gc' in cn_data.columns:
            cn_data['gc'].fillna(0.5)
        else:
            cn_data['gc'] = 0.5
            
        metrics_data.set_index('cell_id', inplace=True)

        cn_unstacked = cn_data \
                        .set_index(['chr', 'start', 'end', 'cell_id'])[['gc', 'reads','copy','state']] \
                        .unstack(level='cell_id').transpose()

        chrom, start, end = zip(*cn_unstacked.columns.tolist())
                
        super().__init__(
            X=cn_unstacked.loc['copy'].fillna(0).values,
            obs={'cell_id': cn_unstacked.loc['copy'].index.tolist()},
            var={'chr': chrom, 'start': start, 'end': end},
            layers={
                'copy': cn_unstacked.loc['copy'].fillna(0).values,
                'gc': cn_unstacked.loc['gc'].fillna(0.5).values,
                'reads': cn_unstacked.loc['reads'].fillna(0).values,
                'state': cn_unstacked.loc['state'].values,
            }
        )

        self.obs.set_index('cell_id', inplace=True)
        self.var.set_index(['chr', 'start', 'end'], inplace=True)
        
        # Make sure we are able to merge
        non_int = set(cn_data['cell_id']).symmetric_difference(metrics_data.index.values)
        if len(non_int) != 0:
            print(f"\x1b[1;33m<WARNING: Cell ids in CN do not match cell ids in metrics: {non_int[0:10]}...>\x1b[0m")
        
        self.obs = self.obs.merge( metrics_data, how='left', left_index=True, right_index=True )
        
        self.set_X('copy', inplace=True)
        