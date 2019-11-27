#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse the config file.
"""

import configparser

def parse_config_file(run_id):
    config_file_path = "{0}/config.txt".format(run_id)

    config = configparser.ConfigParser()
    config.read(config_file_path)

    #GENERAL
    workflow = config.get('GENERAL', 'workflow')
    #DATABASE_PATHS
    root_db = config.get('DATABASE_PATHS','root_db')
    data_base = config.get('DATABASE_PATHS','data_base')
    mnx_folder = config.get('DATABASE_PATHS', 'mnx_folder')
    mnx_rxn = config.get('DATABASE_PATHS','mnx_rxn')
    mnx_cpd = config.get('DATABASE_PATHS', 'mnx_cpd')
    mnx_cpd_prop = config.get('DATABASE_PATHS', 'mnx_cpd_prop')
   
    #PATHS_IN_RUN
    base = config.get('PATHS_IN_RUN', 'base')
    root_path = config.get('PATHS_IN_RUN', 'root_path')
    networks_folder = config.get('PATHS_IN_RUN', 'networks_folder')
    annotation_output_folder = config.get('PATHS_IN_RUN',
                                          'annotation_output_folder')
    orthology_output_folder = config.get('PATHS_IN_RUN',
                                         'orthology_output_folder')
    external_folder = config.get('PATHS_IN_RUN', 'external_folder')
    orthology_model_folder = config.get('PATHS_IN_RUN',
                                        'orthology_model_folder')
    orthofinder_workdir = config.get('PATHS_IN_RUN', 'orthofinder_workdir')
    orthofinder_output = config.get('PATHS_IN_RUN', 'orthofinder_output')
    #dict_gene = config.get('PATHS_IN_RUN', 'dict_gene')
    annotation_folder = config.get('PATHS_IN_RUN', 'annotation_folder')
    curation_data_folder = config.get('PATHS_IN_RUN', 'curation_data_folder')
    genomic_folder = config.get('PATHS_IN_RUN', 'genomic_folder')
    wiki_pages = config.get('PATHS_IN_RUN', 'wiki_pages')
    report_dir = config.get('PATHS_IN_RUN', 'report_dir')
    askomics = config.get('PATHS_IN_RUN', 'askomics')
    faa_study = config.get('PATHS_IN_RUN', 'faa_study')
    gbk_study = config.get('PATHS_IN_RUN', 'gbk_study')
    artefacts = config.get('PATHS_IN_RUN', 'artefacts')
    pathwaytools_output = config.get('PATHS_IN_RUN', 'pathwaytools_output')
    seeds = config.get('PATHS_IN_RUN', 'seeds')
    seeds_artefacts = config.get('PATHS_IN_RUN', 'seeds_artefacts')
    targets = config.get('PATHS_IN_RUN', 'targets')
    meneco_seeds = config.get('PATHS_IN_RUN', 'meneco_seeds')
    meneco_original_output = config.get('PATHS_IN_RUN',
                                        'meneco_original_output')
    meneco_solution = config.get('PATHS_IN_RUN', 'meneco_solution')
    draft = config.get('PATHS_IN_RUN', 'draft')

    #TOOL_PATHS
    programs = config.get('TOOL_PATHS', 'programs')
    padmet_u = config.get('TOOL_PATHS', 'padmet_u')
    reaction_to_add_delete = config.get('TOOL_PATHS', 'reaction_to_add_delete')
    new_reaction_data = config.get('TOOL_PATHS', 'new_reaction_data')

    #TOOL_PARAMETERS
    orthology_method = config.get('TOOL_PARAMETERS', 'orthology_method')
    annotation_method = config.get('TOOL_PARAMETERS', 'annotation_method')
    gap_filling_method = config.get('TOOL_PARAMETERS', 'gap_filling_method')
    all_rxn = config.get('TOOL_PARAMETERS', 'all_rxn')
    with_artefacts = config.get('TOOL_PARAMETERS', 'with_artefacts')
    pwytools_installed = config.get('TOOL_PARAMETERS', 'pwytools_installed')
    no_orphan = config.get('TOOL_PARAMETERS', 'no_orphan')
    cutoff = config.get('TOOL_PARAMETERS', 'cutoff')
    remove_ortho_workdir = config.get('TOOL_PARAMETERS', 'remove_ortho_workdir')
    to_map = config.get('TOOL_PARAMETERS', 'to_map')
    lvl = config.get('TOOL_PARAMETERS', 'lvl')
    
    config_data = {'workflow': workflow, 'data_base': data_base,
                   'mnx_folder': mnx_folder, 'mnx_rxn': mnx_rxn,
                   'mnx_cpd': mnx_cpd, 'mnx_cpd_prop': mnx_cpd_prop,
                   'base': base, 'networks_folder': networks_folder,
                   'annotation_output_folder': annotation_output_folder,
                   'orthology_output_folder': orthology_output_folder,
                   'external_folder': external_folder,
                   'orthology_model_folder': orthology_model_folder,
                   'annotation_folder': annotation_folder,
                   'curation_data_folder': curation_data_folder,
                   'genomic_folder': genomic_folder, 'wiki_pages': wiki_pages,
                   'askomics': askomics, 'faa_study': faa_study,
                   'gbk_study': gbk_study, 'artefacts': artefacts,
                   'pathwaytools_output': pathwaytools_output, 'seeds': seeds,
                   'seeds_artefacts': seeds_artefacts, 'targets': targets,
                   'draft': draft, 'programs': programs, 'padmet_u': padmet_u,
                   'orthofinder_workdir': orthofinder_workdir,
                   'orthofinder_output': orthofinder_output,
                   'meneco_seeds': meneco_seeds,
                   'meneco_original_output': meneco_original_output,
                   'meneco_solution': meneco_solution,
                   'reaction_to_add_delete': reaction_to_add_delete,
                   'new_reaction_data': new_reaction_data,
                   'orthology_method': orthology_method,
                   'annotation_method': annotation_method,
                   'gap_filling_method': gap_filling_method,
                   'all_rxn': all_rxn, 'with_artefacts': with_artefacts,
                   'pwytools_installed': pwytools_installed,
                   'no_orphan': no_orphan, 'cutoff': cutoff, 
                   'remove_ortho_workdir': remove_ortho_workdir,
                   'to_map': to_map, 'lvl': lvl}

    return config_data
