import importlib.util
import os
import time
import functools
import json
from django.http import JsonResponse, HttpRequest, HttpResponse
from django.views.decorators.csrf import csrf_exempt
import pandas as pd
import numpy as np
from collections import Counter
import scanpy as sc
import anndata as an
from . import API
from lifelines.datasets import load_waltons
from collections import Counter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
import math
#####################################################################################################################################
#  functions for specific requests
#####################################################################################################################################

#####################################################################################################################################
#  request handling function
#####################################################################################################################################



@csrf_exempt
def index_of_proportion(request: HttpRequest,
        input_data_all : str):
    """
    response = {
        "propotion":{
            array
        }
    }
    """
    c = API.DatabaseAPI("bulk")
    my_dict_c = c.query_collection_obs()
    my_df_c = pd.DataFrame(my_dict_c)

    input_list_all = input_data_all.split("&")
    input_cell = input_list_all[0].split("=")[1]
    if ("," in input_cell):
        input_cell_list = input_cell.split(",")
    else:
        input_cell_list = [input_cell]
    input_data = input_list_all[1].split("=")[1]
    if ("," in input_data):
        input_list = input_data.split(",")
    else:
        input_list = [input_data]
    ref = input_list_all[3].split("=")[1]
    

    if len(input_list) > 8 :
        response = {"error" : "Too many datasets. You can select no more than eight datasets."}
        return JsonResponse(response)

    my_df_d = my_df_c[my_df_c["select"].isin(input_list)]
    select = my_df_d["select"].tolist()
    if ref == "EPIC":
        columns_list = ["EPIC_cellFractions."+i for i in input_cell_list]
    elif ref == "LM":
        columns_list = ["LM_"+i for i in input_cell_list]
    elif ref == "QS":
        columns_list = ["QS_"+i for i in input_cell_list]
    else:
        response = {"error" : "reference error"}
        return JsonResponse(response)

    my_df_d = my_df_d.loc[:,columns_list]
    if input_list_all[2].split("=")[1] == "true":
        my_df_d = my_df_d.div(my_df_d.sum(axis = 1), axis='rows')
        my_df_d = my_df_d.fillna(0)
    my_df_d = np.round(my_df_d, 4)
    my_df_d.insert(0,"labels",select)
    my_df_array =  np.array(my_df_d.T)


    response = {
        "propotion": my_df_array.tolist(),
    }
    return JsonResponse(response)


def index_of_expression(request: HttpRequest,
        summary : str):
    """
    response = {
        "expression":{
            array
        }
        "gene" : {
            list
        }
    }
    """
    input_data_all = summary.split("&")[0] + "&" + summary.split("&")[1] + "&norm=false"
    c = API.DatabaseAPI("bulk")
    my_dict_c = c.query_collection_obs()
    my_df_c = pd.DataFrame(my_dict_c)

    input_list_all = input_data_all.split("&")
    input_cell = input_list_all[0].split("=")[1]
    if ("," in input_cell):
        input_cell_list = input_cell.split(",")
    else:
        input_cell_list = [input_cell]
    input_data = input_list_all[1].split("=")[1]
    input_list = input_data.split(",")
    gene_list = summary.split("&")[2].split("=")[1]
    if ("," in gene_list):
        gene_list = gene_list.split(",")
    else:
        gene_list = [gene_list]

    if input_list_all[0] == "type=":
        expression = c.query_collection_gene_X_var_by_obs(gene_list)
        expression = pd.DataFrame(expression)
        my_dict_c = c.query_collection_obs()
        my_df_c = pd.DataFrame(my_dict_c)
        select_part = my_df_c.loc[my_df_c["select"].isin(input_list),:]
        expression.columns = my_df_c["cellID"].tolist()
        expression_part = expression.loc[:,select_part["cellID"].tolist()]
        expression_array = np.round(np.log(np.array(expression_part)),4)
        expression_df = pd.DataFrame(expression_array).T
        expression_df.insert(0,"select",select_part["select"].tolist())
        expression_array =  np.array(expression_df.T)
        response = {
            "expression": expression_array.tolist(),
            "genes":gene_list
        }
        return JsonResponse(response)




    ref = summary.split("&")[3].split("=")[1]


    if len(gene_list) > 5 :
        response = {"error" : "Too many genes. You can select no more than five genes."}
        return JsonResponse(response)

    if len(input_list) > 8 :
        response = {"error" : "Too many datasets. You can select no more than eight datasets."}
        return JsonResponse(response)



    my_df_d = my_df_c[my_df_c["select"].isin(input_list)]

    select = my_df_d["select"].tolist()

    if ref == "EPIC":
        columns_list = ["EPIC_cellFractions."+i for i in input_cell_list]
        ref_database = API.DatabaseAPI("ref")
    elif ref == "LM":
        columns_list = ["LM_"+i for i in input_cell_list]
        ref_database = API.DatabaseAPI("LM_ref")
    elif ref == "QS":
        columns_list = ["QS_"+i for i in input_cell_list]
        ref_database = API.DatabaseAPI("QS_ref")
    else:
        response = {"error" : "reference error"}
        return JsonResponse(response)

    cellID = my_df_d["cellID"].tolist()
    my_df_d = my_df_d.loc[:,columns_list]
    if input_list_all[2].split("=")[1] == "true":
        my_df_d = my_df_d.div(my_df_d.sum(axis = 1), axis='rows')
    my_df_d.index = cellID
    my_df_t = my_df_d.T

    
    genes = ref_database.query_collection_var()["geneSymbol"]
    gene_list = [i for i in gene_list if i in genes]
    
    gg = ref_database.query_collection_gene_X_var_by_obs(gene_list)
    gg = pd.DataFrame(gg)


    gg.columns = ref_database.query_collection_obs()["celltype"]
    gg_mean = pd.DataFrame(gg.T.mean(axis = 1))
    gg_mean = gg_mean.loc[input_cell_list,:]
    expression = my_df_t.multiply(gg_mean.values)
    expression_t = expression.T
    expression_t = pd.DataFrame(np.log1p(np.array(expression_t)))
    expression_t = np.round(expression_t, 4)
    expression_t.insert(0,"select",select)
    expression_array =  np.array(expression_t.T)

    response = {
        "expression": expression_array.tolist(),
        "genes":gene_list
    }
    return JsonResponse(response)


def index_of_survival(request: HttpRequest,
        all_parameter : str):
    """
    response = {
        data = 
    }
    """ 
    mm = all_parameter.split("&")

    st = mm[0].split("=")[1]
    if ("," in mm[1].split("=")[1]):
        ct = mm[1].split("=")[1].split(",")
    else:
        ct = [mm[1].split("=")[1]]

    b = API.DatabaseAPI("tcga")
    my_dict_b = b.query_collection_obs()
    my_df_b = pd.DataFrame(my_dict_b)
    select_part = my_df_b.loc[my_df_b["primary_disease"].isin(ct),:]

    if len(ct) > 8 :
        response = {"error" : "Too many datasets. You can select no more than eight datasets."}
        return JsonResponse(response)

    ref = mm[5].split("=")[1]
    if ("," in mm[2].split("=")[1]):
        cell = mm[2].split("=")[1].split(",")
    else:
        cell = [mm[2].split("=")[1]]
    up = mm[3].split("=")[1]
    dn = mm[4].split("=")[1]

    select = select_part["primary_disease"].tolist()
    if ref == "EPIC":
        columns_list = ["EPIC_cellFractions."+i for i in cell]
        ref = API.DatabaseAPI("ref")
    elif ref == "LM":
        columns_list = ["LM_"+i for i in cell]
        ref = API.DatabaseAPI("LM_ref")
    elif ref == "QS":
        columns_list = ["QS_"+i for i in cell]
        ref = API.DatabaseAPI("QS_ref")
    else:
        response = {"error" : "reference error"}
        return JsonResponse(response)
    cellID = select_part["cellID"].tolist()
    my_df_d = select_part.loc[:,columns_list]
    my_df_d.index = cellID
    my_df_t = my_df_d.T

        
    genes = ref.query_collection_var()["geneSymbol"]
    
    gg = ref.query_collection_gene_X_var_by_obs(genes)
    gg = pd.DataFrame(gg)

    gg.columns = ref.query_collection_obs()["celltype"]
    gg_mean = pd.DataFrame(gg.T.mean(axis = 1))
    gg_mean = gg_mean.loc[cell,:]
    expression = my_df_t.multiply(gg_mean.values)
    expression_t = expression.T
    expression_t = pd.DataFrame(expression_t.sum(axis = 1),columns = ["sum"])
    expression_t = expression_t.sort_values(by=["sum"],ascending = False)
    number = expression_t.shape[0]
    number1 = int(number/100 * (100-int(up)))
    number2 = int(number/100 * (100-int(dn)))
    samples = expression_t.index.tolist()
    sample = []
    for each in samples:
        names = each.split(".")
        sample.append(names[0]+"."+names[1]+"."+names[2])

    up_sample = sample[:number1]
    dn_sample = sample[number2+1:]

    matches = {"Dead":1,"Alive":0,"-":0}
    a = API.DatabaseAPI("survival")

    #up part
    my_dict_a = a.query_collection_obs()
    my_df_a = pd.DataFrame(my_dict_a)
    my_df_a = my_df_a.loc[my_df_a["sample"].isin(up_sample),:]
    OSEVENT = my_df_a["OSEVENT"].tolist()
    E = [matches[i] for i in OSEVENT]
    if st == "OS":
        T = my_df_a["OSDAY"].tolist()
    else:
        T = my_df_a["RFSDAY"].tolist()

    E_end_up = [E[i] for i in range(len(T)) if T[i] != "-"]
    T_end = [T[i] for i in range(len(T)) if T[i] != "-"]
    T_end = list(map(float, T_end))
    T_end_up =  list(map(lambda x: round(x/30,2), T_end))
    kmf = KaplanMeierFitter()
    kmf.fit(T_end_up, E_end_up)

    sf=kmf.survival_function_.T
    xa=sf.columns.tolist()
    y1a =  list(map(lambda x: round(x,3), sf.values[0].tolist()))
    ci=kmf.confidence_interval_survival_function_.T.values
    y2a =  list(map(lambda x: round(x,3), ci[1].tolist()))
    y3a =  list(map(lambda x: round(x,3), ci[0].tolist()))
    xca = [T_end_up[i] for i in range(len(T_end_up)) if E_end_up[i]==0]
    xca = list(map(float, xca))
    yca =  list(map(lambda x: round(x,3), kmf.survival_function_at_times(xca).tolist()))

    #dn part
    my_dict_a = a.query_collection_obs()
    my_df_a = pd.DataFrame(my_dict_a)
    my_df_a = my_df_a.loc[my_df_a["sample"].isin(dn_sample),:]

    OSEVENT = my_df_a["OSEVENT"].tolist()
    E = [matches[i] for i in OSEVENT]
    if st == "OS":
        T = my_df_a["OSDAY"].tolist()
    else:
        T = my_df_a["RFSDAY"].tolist()

    E_end_dn = [E[i] for i in range(len(T)) if T[i] != "-"]
    T_end = [T[i] for i in range(len(T)) if T[i] != "-"]
    T_end = list(map(float, T_end))
    T_end_dn =  list(map(lambda x: round(x/30,2), T_end))
    kmf = KaplanMeierFitter()
    kmf.fit(T_end_dn, E_end_dn)

    sf=kmf.survival_function_.T
    xb=sf.columns.tolist()
    y1b =  list(map(lambda x: round(x,3), sf.values[0].tolist()))
    ci=kmf.confidence_interval_survival_function_.T.values
    y2b =  list(map(lambda x: round(x,3), ci[1].tolist()))
    y3b =  list(map(lambda x: round(x,3), ci[0].tolist()))
    xcb=[T_end_dn[i] for i in range(len(T_end_dn)) if E_end_dn[i]==0]
    xcb = list(map(float, xcb))
    ycb =  list(map(lambda x: round(x,3), kmf.survival_function_at_times(xcb).tolist()))

    results = logrank_test(T_end_up, T_end_dn, event_observed_A=E_end_up, event_observed_B=E_end_dn)
    pValues1 = float(results.summary["p"].values)

    dfA = pd.DataFrame({'E': E_end_up, 'T': T_end_up, 'groupA': 1})
    dfB = pd.DataFrame({'E': E_end_dn, 'T': T_end_dn, 'groupA': 0})
    df = pd.concat([dfA, dfB])
    cph = CoxPHFitter().fit(df, 'T', 'E')
    pValues2 = float(cph.summary["p"].values)

    response = {"data":
            [
        {
            "pValues1" : pValues1,
            "pValues2" : pValues2
        },
        {
            "line": {
                "dash": "solid",
                "color": "red",
                "shape": "hv",
                "width": 2
            },
            "mode": "lines",
            "name": "",
            "type": "scatter",
            "x": xa,
            "y": y1a,
            "xaxis": "x1",
            "yaxis": "y1",
            "showlegend": False
        },
        {
            "line": {
                "dash": "dash",
                "color": "red",
                "shape": "hv",
                "width": 2
            },
            "mode": "lines",
            "name": "",
            "type": "scatter",
            "x": xa,
            "y": y2a,
            "xaxis": "x1",
            "yaxis": "y1",
            "showlegend": False
        },
        {
            "line": {
                "dash": "dash",
                "color": "red",
                "shape": "hv",
                "width": 2
            },
            "mode": "lines",
            "name": "",
            "type": "scatter",
            "x": xa,
            "y": y3a,
            "xaxis": "x1",
            "yaxis": "y1",
            "showlegend": False
        },
        {
            "mode": "markers",
            "name": "",
            "text": "",
            "type": "scatter",
            "x": xca,
            "y": yca,
            "xaxis": "x1",
            "yaxis": "y1",
            "marker": {
                "size": 10,
                "color": "black",
                "symbol": "cross-thin-open",
                "opacity": 1,
                "sizeref": 1,
                "sizemode": "area"
            },
            "showlegend": False
        },
        { 
            "line" : { 
                "dash": "solid", 
                "color": "blue", 
                "shape": "hv", 
                "width": 2 
            }, 
            "mode": "lines", 
            "name": "", "type": 
            "scatter", 
            "x":xb, 
            "y": y1b, 
            "xaxis": "x1", 
            "yaxis": "y1", 
            "showlegend": False
        }, 
        { 
            "line": { 
                "dash": "dash", 
                "color": "blue", 
                "shape": "hv", 
                "width": 2 
            }, 
            "mode": "lines",
            "name": "", 
            "type": "scatter", 
            "x": xb, 
            "y": y2b, 
            "xaxis": "x1", 
            "yaxis": "y1", 
            "showlegend": False 
        }, 
        { 
            "line": { 
                "dash": "dash", 
                "color": "blue", 
                "shape": "hv", 
                "width": 2 
            }, 
            "mode": "lines", 
            "name": "", 
            "type": "scatter", 
            "x": xb, 
            "y": y3b, 
            "xaxis": "x1", 
            "yaxis": "y1", 
            "showlegend": False
        }, 
        { 
            "mode": "markers", 
            "name": "", 
            "text": "", 
            "type": "scatter", 
            "x": xcb, 
            "y": ycb, 
            "xaxis": "x1", 
            "yaxis": "y1", 
            "marker": { 
                "size": 10, 
                "color": "black", 
                "symbol": "cross-thin-open", 
                "opacity": 1, 
                "sizeref": 1, 
                "sizemode": "area" 
                }, 
            "showlegend": False
        }
    ]
    }
    return JsonResponse(response)


