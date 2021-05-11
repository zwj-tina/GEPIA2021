import os
from typing import Optional, Union, Dict, Any, List

import numpy as np
import pandas as pd
from pymongo import MongoClient, ASCENDING

# the API to write and read
class DatabaseAPI(object):

    def __init__(
            self,
            target_db: Optional[str] = None
    ):
        self.target_db = target_db
        self.client = MongoClient()
        self.db = self.client[self.target_db]

    # =======================
    #   write methods
    # =======================
    # write the obs_by_var into database
    def write_collection_matrix_obs_by_var(
            self,
            obs_by_var_matrix: Optional[Union[np.ndarray]] = None,
            overwrite: bool = False,
            verbose: bool = True
    ):

        collection_name = 'expression_matrix_obs_by_var'
        index_key = 'obs_id'
        # if exist, choose overwrite or not
        if collection_name in self.db.list_collection_names() and not overwrite:
            return_message = "collection %s exists, set overwrite=True if you want to overwrite" % collection_name
            if verbose:
                print(return_message)
            return return_message
        if collection_name in self.db.list_collection_names() and overwrite:
            self.db.drop_collection(collection_name)

        documents = []
        num_of_document = 0
        for i in range(obs_by_var_matrix.shape[0]):
            document = {
                "obs_id": i,
                "obs_value": obs_by_var_matrix[i].tolist()
            }
            documents.append(document)
            num_of_document += 1

        collection = self.db[collection_name]
        collection.insert_many(documents, bypass_document_validation=True)

        # create index to avoid document traversal in query.
        collection.create_index([(index_key, ASCENDING)])

        return_message = 'successfully loaded in %d document in collection %s and created index on field %s ' % (
            num_of_document, collection_name, index_key)
        return return_message
    
    def write_collection_gene_matrix_var_by_obs(
            self,
            var_by_obs_matrix: Optional[Union[np.ndarray]] = None,
            genes = [],
            overwrite: bool = False,
            verbose: bool = True
    ):

        collection_name = 'expression_matrix_gene_var_by_obs'
        index_key = 'var_id'
        # if exist, choose overwrite or not
        if collection_name in self.db.list_collection_names() and not overwrite:
            return_message = "collection %s exists, set overwrite=True if you want to overwrite" % collection_name
            if verbose:
                print(return_message)
            return return_message
        if collection_name in self.db.list_collection_names() and overwrite:
            self.db.drop_collection(collection_name)

        documents = []
        num_of_document = 0
        for i in range(var_by_obs_matrix.shape[0]):
            document = {
                "var_id": genes[i],
                "var_value": var_by_obs_matrix[i].tolist()
            }
            documents.append(document)
            num_of_document += 1

        collection = self.db[collection_name]
        collection.insert_many(documents, bypass_document_validation=True)

        # create index to avoid document traversal in query.
        collection.create_index([(index_key, ASCENDING)])

        return_message = 'successfully loaded in %d document in collection %s and created index on field %s ' % (
            num_of_document, collection_name, index_key)
        return return_message
    # write var into database
    def write_collection_var(
            self,
            var: Optional[Dict[str, list]] = None,
            overwrite: bool = False,
            verbose: bool = True
    ):
        collection_name = 'var'
        # if exist, choose overwrite or not
        if collection_name in self.db.list_collection_names() and not overwrite:
            return_message = "collection %s exists, set overwrite=True if you want to overwrite" % collection_name
            if verbose:
                print(return_message)
            return return_message
        if collection_name in self.db.list_collection_names() and overwrite:
            self.db.drop_collection(collection_name)

        documents = []

        for key in var.keys():
            document = {
                'key': key,
                'value': var[key]
            }
            documents.append(document)

        collection = self.db[collection_name]
        collection.insert_many(documents)

        # no need of index

        return_message = 'successfully loaded in collection %s ' % collection_name
        return return_message

    # write obs into database
    def write_collection_obs(
            self,
            obs: Optional[Dict[str, list]] = None,
            overwrite: bool = False,
            verbose: bool = True
    ):

        collection_name = 'obs'
        # if exist, choose overwrite or not
        if collection_name in self.db.list_collection_names() and not overwrite:
            return_message = "collection %s exists, set overwrite=True if you want to overwrite" % collection_name
            if verbose:
                print(return_message)
            return return_message
        if collection_name in self.db.list_collection_names() and overwrite:
            self.db.drop_collection(collection_name)

        documents = []

        for key in obs.keys():
            document = {
                'key': key,
                'value': obs[key]
            }
            documents.append(document)

        collection = self.db[collection_name]
        collection.insert_many(documents)

        return_message = 'successfully loaded in collection %s ' % collection_name
        return return_message
    # =======================
    #   read methods
    # =======================
    # read the obs_by_var from the database,parameter:the index of obs  type int
    def read_collection_X_obs_by_var(
            self,
            obs_id: Optional[int] = None
    ) -> list:

        """
        returns the obs_value(list) of a given obs_id
        """
        collection_name = 'expression_matrix_obs_by_var'
        collection = self.db[collection_name]

        # use find_one() but not find()
        obs = collection.find_one({'obs_id': obs_id})
        obs_value = obs['obs_value']

        return obs_value
    
    def read_collection_gene_X_var_by_obs(
            self,
            gene_id = None,
    ) -> list:

        """
        returns the obs_value(list) of a given obs_id
        """
        collection_name = 'expression_matrix_gene_var_by_obs'
        collection = self.db[collection_name]

        # use find_one() but not find()
        var = collection.find_one({'var_id': gene_id})
        var_value = var['var_value']

        return var_value

 
    # read the var from the database,parameter:key word  type:str
    def read_collection_var(
            self,
            key: Optional[str] = None
    ) -> list:

        """
        returns the value(list) of a given key in var
        """
        collection_name = 'var'
        collection = self.db[collection_name]

        v = collection.find_one({'key': key})
        value = v['value']

        return value

    # read the obs from the database,parameter:key word  type:str
    def read_collection_obs(
            self,
            key: Optional[str] = None
    ) -> list:

        """
        returns the value(list) of a given key in obs
        """
        collection_name = 'obs'
        collection = self.db[collection_name]

        v = collection.find_one({'key': key})
        value = v['value']
        return value
    def query_collection_X_obs_by_var(
            self,
            obs_ids: Optional[List[int]] = None
    ) -> np.ndarray:

        """
        returns a numpy array
        """
        obs_values = []
        for obs_id in obs_ids:
            obs_value = self.read_collection_X_obs_by_var(obs_id)
            obs_values.append(obs_value)

        obs_values = np.array(obs_values)
        return obs_values
    
    def query_collection_gene_X_var_by_obs(
            self,
            gene_ids = []
    ) -> np.ndarray:

        """
        returns a numpy array
        """
        var_values = []
        for var_id in gene_ids:
            var_value = self.read_collection_gene_X_var_by_obs(var_id)
            var_values.append(var_value)

        var_values = np.array(var_values)
        return var_values

    # query the var_by_obs from the database,parameter:the index of var type list
   
    def query_collection_obs(
            self
    ) -> [Dict[str, Any]]:

        collection_name = 'obs'
        collection = self.db[collection_name]

        values = dict()
        cursor = collection.find()
        for each in cursor:
            name = each['key']
            values[name] = each['value']
        return values

    # query the all var from the database,no parameter
    def query_collection_var(
            self
    ) -> [Dict[str, Any]]:

        collection_name = 'var'
        collection = self.db[collection_name]

        values = dict()
        cursor = collection.find()
        for each in cursor:
            name = each['key']
            values[name] = each['value']
        return values
    def get_collection_X_obs_by_var(
            self
    ) -> np.ndarray:

        """
        returns the matrix in np.array format
        """
        collection_name = 'expression_matrix_obs_by_var'
        collection = self.db[collection_name]

        obs_by_var_matrix = None

        obs_ids = collection.count()
        obs_by_var_matrix = self.query_collection_X_obs_by_var(range(obs_ids))

        return obs_by_var_matrix
    
    def get_collection_gene_X_var_by_obs(
            self
    ) -> np.ndarray:

        """
        returns the matrix in np.array format
        """
        collection_name = 'expression_matrix_gene_var_by_obs'
        collection = self.db[collection_name]

        var_by_obs_matrix = None
        
        gene_ids = []
        cursor = collection.find()
        for each in cursor:
            gene_ids.append(each['var_id'])
        
        var_by_obs_matrix = self.query_collection_gene_X_var_by_obs(gene_ids)

        return var_by_obs_matrix

