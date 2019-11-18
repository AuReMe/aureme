"""
usage:
    aureme --init=ID [-v]
    aureme <command> [<args>...]
    aureme -R --run=ID
    aureme --version
    aureme --installPWT=PWT_path [--ptools=ptools_path]
    aureme --uninstallPWT

options:
    -h --help     Show help.
    -R     Give access from container.
    --init=ID    Create folder structure.
    --run=ID    Pathname to the comparison workspace.
    -v     Verbose.

The subcommands are:
    check    Check inputs validity
    reconstruction    Run Pathway Tools
    orthology    Run Orthofinder for crossing orthology between species
    draft    Merge all networks (from Pathway Tools and Orthology)

    workflow    Run Check, Pathway Tools, Orthofinder and Merging of all networks
    analysis    Analyze results
    compare    Compare group of species

See 'aureme <command> -h' for more information on a specific command.
"""
import configparser
import csv
import docopt
import eventlet
import mpwt
import os
import os.path
import pkg_resources
import re
import requests
import shutil
import subprocess
import sys

import aureme

release_on_gitlab = "https://raw.githubusercontent.com/AuReMe/aureme/master/release.txt"


def main(args=None):
    args = docopt.docopt(__doc__, options_first=True)

    command = args.pop('<command>')
    command_args = args.pop('<args>')


    if args["--init"]:
        run_id = args["--init"]
        create_run(run_id)
        get_full_right(run_id)
        return

    # Add permission to all folder in run_id, usefull because all cmd exec from container are root based.
    if args['-R']:
        run_id = args["--run"]
        get_full_right(run_id)
        return

    if args['--installPWT']:
        installing_pwt(args['--installPWT'], args['--ptools'])
        return

    if args['--uninstallPWT']:
        uninstalling_pwt()
        return

    if args["--version"]:
        online_version = get_version()
        current_version = pkg_resources.get_distribution("metage2metabo").version
        if online_version:
            print("You are using the version %s, the latest is %s" %(current_version, online_version))
        else:
            print('No internet connection. Skip checking aureme version.')
        return

    if command:
        if command not in ['workflow', 'check', 'reconstruction', 'orthology', 'draft', 'analysis', 'compare']:
            sys.exit(command + ' not a valid command: workflow, check, reconstruction, orthology, draft, analysis, compare.')

        if '-h' in command_args:
            getattr(aureme, command).command_help()
            sys.exit()

        # Add command to command_args to be parse by docopt.
        command_args.insert(0,command)

        if command == 'workflow':
            aureme.workflow.workflow_parse_args(command_args)

        elif command == 'check':
            aureme.check.check_parse_args(command_args)

        elif command == 'reconstruction':
            aureme.reconstruction.reconstruction_parse_args(command_args)
        
        elif command == 'orthology':
            aureme.orthology.orthology_parse_args(command_args)

        elif command == 'draft':
            aureme.draft.draft_parse_args(command_args)

        elif command == 'analysis':
            aureme.analysis.analysis_parse_args(command_args)

        elif command == 'compare':
            aureme.compare.compare_parse_args(command_args)


def create_run(run_id):
    """
    Create the folders of Aureme workspace.
    """
    if os.path.isdir('{0}'.format(run_id)):
        print("Run '%s' already exist, remove this folder manually before"
              %run_id)
    else:
        print('------>RUNNING STEP : Working directory initialization: run %s'
              %run_id)
        os.mkdir('{0}'.format(run_id))
        all_folders = ['analysis', 'analysis/flux_analysis', 'analysis/report',
                       'analysis/topological_analysis', 'analysis/askomics',
                       'analysis/wiki_pages',
                       'annotation_based_reconstruction',
                       'database',
                       'gap_filling', 'gap_filling/original_output',
                       'genomic_data',
                       'growth_medium',
                       'manual_curation', 'manual_curation/template',
                       'networks', 'networks/external_network/',
                       'networks/output_annotation_based_reconstruction',
                       'networks/output_annotation_based_reconstruction/pathwaytools',
                       'networks/output_orthology_based_reconstruction',
                       'networks/output_orthology_based_reconstruction/orthofinder',
                       'orthology_based_reconstruction',
                       'orthology_based_reconstruction/orthofinder_wd',
                       'targets_compounds']
        for folder in all_folders:
            print('creating folder {0}/{1}'.format(run_id, folder))
            os.mkdir("{0}/{1}".format(run_id, folder))

        create_default_file(run_id)
        
        config_file_path = '{0}/config.txt'.format(run_id)
        create_config_file(config_file_path, run_id)


def create_config_file(config_file_path, run_id):
    config = configparser.RawConfigParser()
    config.add_section('GENERAL')
    config.set('GENERAL', 'workflow', 'AuReMe')
    config.set('GENERAL', 'run_id', os.path.basename(run_id))
    #config.set('GENERAL', 'new_network', '%(network)s')
    
    config.add_section('DATABASE_PATHS')
    config.set('DATABASE_PATHS', '#data_base', '%(run_id)s/database/XXX')
    config.set('DATABASE_PATHS', 'data_base',
               '/home/data/database/BIOCYC/METACYC/23.0/metacyc_23.0')
    config.set('DATABASE_PATHS', 'mnx_folder', '/home/data/database/MNX/2018')
    config.set('DATABASE_PATHS', 'mnx_rxn', '%(mnx_folder)s/reac_xref.tsv')
    config.set('DATABASE_PATHS', 'mnx_cpd', '%(mnx_folder)s/chem_xref.tsv')
    config.set('DATABASE_PATHS', 'mnx_cpd_prop', '%(mnx_folder)s/chem_prop.tsv')
   
    config.add_section('PATHS_IN_RUN')
    config.set('PATHS_IN_RUN', 'base', run_id)
    config.set('PATHS_IN_RUN', 'networks_folder', '%(base)s/networks')
    config.set('PATHS_IN_RUN', 'annotation_output_folder',
               '%(base)s/networks/output_annotation_based_reconstruction')
    config.set('PATHS_IN_RUN', 'orthology_output_folder',
               '%(base)s/networks/output_orthology_based_reconstruction')
    config.set('PATHS_IN_RUN', 'external_folder',
               '%(base)s/networks/external_network')
    config.set('PATHS_IN_RUN', 'orthology_model_folder',
               '%(base)s/orthology_based_reconstruction')
    config.set('PATHS_IN_RUN', 'dict_gene', '%(model_folder)s/dict_genes.txt')
    config.set('PATHS_IN_RUN', 'annotation_folder',
               '%(base)s/annotation_based_reconstruction')
    config.set('PATHS_IN_RUN', 'curation_data_folder',
               '%(base)s/manual_curation')
    config.set('PATHS_IN_RUN', 'genomic_folder', '%(base)s/genomic_data')
    config.set('PATHS_IN_RUN', 'wiki_pages', '%(base)s/analysis/wiki_pages')
    config.set('PATHS_IN_RUN', 'report_dir', '%(base)s/analysis/report')
    config.set('PATHS_IN_RUN', 'askomics', '%(base)s/analysis/askomics')
    config.set('PATHS_IN_RUN', 'faa_study', '%(genomic_folder/run_id)s.faa')
    config.set('PATHS_IN_RUN', 'gbk_study', '%(genomic_folder/run_id)s.gbk')
    config.set('PATHS_IN_RUN', 'artefacts', 'growth_medium/artefacts')
    config.set('PATHS_IN_RUN', 'pathwaytools_output',
               '%(networks_folder)s/output_annotation_based_reconstruction/pathwaytools/output_pathwaytools')
    config.set('PATHS_IN_RUN', 'seeds','growth_medium/seeds')
    config.set('PATHS_IN_RUN', 'seeds_artefacts',
               'growth_medium/seeds_artefacts')
    config.set('PATHS_IN_RUN', 'targets', 'targets_compounds/targets')
    config.set('PATHS_IN_RUN', 'draft', '%(networks_folder)s/draft')    
    
    config.add_section('TOOL_PATHS')
    config.set('TOOL_PATHS', 'programs', '/programs')
    config.set('TOOL_PATHS', 'padmet_u', '%(programs)s/padmet-utils')
    config.set('TOOL_PATHS', 'orthofinder_workdir',
               '%(orthology_model_folder)s/orthofinder_wd')
    config.set('TOOL_PATHS', 'orthofinder_output',
               '%(orthofinder_workdir)s/Orthologues')
    config.set('TOOL_PATHS', 'meneco_seeds', '%(seeds)s')
    config.set('TOOL_PATHS', 'meneco_original_output',
               '%(base)s/gapfilling/original_output/meneco_output_%(run_is)s.txt')
    config.set('TOOL_PATHS', 'meneco_solution',
               '%(base)s/gapfilling/gapfilling_solution_with_meneco_%(run_id)s.csv')
    config.set('TOOL_PATHS', 'reaction_to_add_delete', 'manual_curation/data/reaction_to_add_delete.csv')
    config.set('TOOL_PATHS', 'new_reaction_data', 'manual_curation/data/reaction_creator.csv')

    config.add_section('TOOL_PARAMETERS')
    config.set('TOOL_PARAMETERS', 'orthology_method', 'orthofinder')
    config.set('TOOL_PARAMETERS', 'annotation_method', 'pathwaytools')
    config.set('TOOL_PARAMETERS', 'gap_filling_method', 'meneco')
    config.set('TOOL_PARAMETERS', 'all_rxn', '-f')
    config.set('TOOL_PARAMETERS', 'with_artefacts', 'TRUE')
    config.set('TOOL_PARAMETERS', 'pwytools_installed', 'FALSE')
    config.set('TOOL_PARAMETERS', 'no_orphan', '--no-orphan')
    config.set('TOOL_PARAMETERS', 'cutoff','0.70')
    config.set('TOOL_PARAMETERS', 'remove_ortho_workdir', 'TRUE')
    config.set('TOOL_PARAMETERS', 'to_map', 'reaction')
    config.set('TOOL_PARAMETERS', 'lvl', '3')
    
    # Writing our configuration file to 'example.cfg'
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)

        
def create_default_file(run_id):
    '''
    Create the default files of Aureme: full_log.txt, log.txt, 
    default_artefacts_metacyc.txt, reaction_creator.csv, and
    reaction_to_add_delete.csv
    '''
    # full_log.txt file creation
    with open('{0}/full_log.txt'.format(run_id), 'w') as full_log:
        full = csv.writer(full_log, delimiter='\t')
        full.writerow(['### FULL LOG ###'])
    # log.txt file creation
    with open('{0}/log.txt'.format(run_id), 'w') as txt_log:
        log = csv.writer(txt_log, delimiter='\t')
        log.writerow(['### LOG ###'])

    # default_artefacts_metacyc.txt file creation
    with open('{0}/{1}/default_artefacts_metacyc.txt'.format(run_id, 'growth_medium'), 'w') as default_artefact:
        artefact = csv.writer(default_artefact, delimiter='\t')
        artefact.writerow(['NADP', 'c'])
        artefact.writerow(['ADP', 'c'])

    # reaction_to_add_delete.csv file creation
    with open('{0}/{1}/{2}/reaction_to_add_delete.csv'.format(run_id, 'manual_curation', 'template'), 'w') as reaction_add:
        react_add = csv.writer(reaction_add, delimiter='\t')
        react_add.writerow(['idRef',	'Comment', 'Action', 'Genes'])
        react_add.writerow(['rxn_id_1', 'Reaction deleted for x reason',
                            'delete'])
        react_add.writerow(['rxn_id_2',	'Reaction added for x reason', 'add',
                            '(gene1 and gene2)'])
        react_add.writerow(['rxn_id_3', 'Reaction added for x reason', 'add'])

    # reaction_creator.csv file creation
    with open('{0}/{1}/{2}/reaction_creator.csv'.format(run_id, 'manual_curation', 'template'), 'w') as reaction_creator:
        react_cr = csv.writer(reaction_creator, delimiter='\t')
        react_cr.writerow(['reaction_id', 'my_rxn'])
        react_cr.writerow(['comment', 'reaction added for X reason'])
        react_cr.writerow(['reversible', 'false'])
        react_cr.writerow(['linked_gene', '(gene_a or gene_b) and gene_c'])
        react_cr.writerow(['#reactant/product',
                           '#stoichio:compound_id:compart'])
        react_cr.writerow(['reactant', '1.0:compound_a:c'])
        react_cr.writerow(['reactant', '2.0:compound_b:c'])
        react_cr.writerow(['product', '1.0:compound_c:c'])
        react_cr.writerow([])
        react_cr.writerow(['reaction_id', 'my_rxn_2'])
        react_cr.writerow(['comment', 'reaction added for X reason'])
        react_cr.writerow(['reversible', 'true'])
        react_cr.writerow(['linked_gene', ''])
        react_cr.writerow(['#reactant/product',
                           '#stoichio:compound_id:compart'])
        react_cr.writerow(['reactant', '1.0:compound_a:c'])
        react_cr.writerow(['reactant', '2.0:compound_d:c'])
        react_cr.writerow(['product', '1.0:compound_c:c'])
        return
            
                  
def get_full_right(name):
    '''
    Get full rigths to the name (file or directory).
    '''
    chmod_cmds = ["chmod", "-R", "777", name]
    subprocess.call(chmod_cmds)

    
def get_version():
    '''
    Get version from Gitlab.
    Check internet connection using requests and eventlet timeout.
    '''
    reg_version = r'^\#+VERSION:([0-9.]*)#+'
    with eventlet.Timeout(2):
        try:
            response = requests.get(release_on_gitlab)
            first_line = response.text.split('\n')[0]
            version = re.match(reg_version,first_line).group(1)
        except eventlet.timeout.Timeout:
            version = None

    return version


def installing_pwt(pwt_path, input_ptools_local_path):
    """
    Install silently Pathway-Tools in /programs.
    After running this function you need to source the bashrc.
    """
    ptools_local_path = '/root'
    if input_ptools_local_path:
        if os.path.isdir(input_ptools_local_path):
            ptools_local_path = input_ptools_local_path
        else:
            print(input_ptools_local_path + ' path does not exist, --ptools must be an existing path.')
            return

    cmd_chmods = ['chmod', 'u+x', pwt_path]
    cmd_installs = [pwt_path, '--InstallDir', '/programs/pathway-tools', '--PTOOLS_LOCAL_PATH', ptools_local_path,
                    '--InstallDesktopShortcuts', '0', '--mode', 'unattended']

    print(' '.join(cmd_chmods))
    subprocess.call(cmd_chmods)
    print(' '.join(cmd_installs))
    subprocess.call(cmd_installs)

    # Add Pathway-Tools in the PATH.
    cmd_echo = '''echo 'export PATH="$PATH:/programs/pathway-tools:"' >> ~/.bashrc'''
    print(cmd_echo)
    subprocess.call(cmd_echo, shell=True)

    print('Now you need to source your bash, run:')
    print('source ~/.bashrc')
    return


def uninstalling_pwt():
    """
    Uninstall Pathway-Tools and can delete ptools-local folder.
    """
    def ask_delete_ptools(ptools_path):
        yes_or_no = input('Delete ptools-local folder (y/n)?')
        if yes_or_no == 'y':
            shutil.rmtree(ptools_path)
            print('Uninstallation of Pahtway-Tools and ptools-local done!')
            return
        elif yes_or_no == 'n':
            print('Uninstallation of Pathway-Tools done!.')
            return
        else:
            print('Wrong command')
            ask_delete_ptools(ptools_path)

    ptools_path = mpwt.find_ptools_path()

    cmd_uninstall = ['/programs/pathway-tools/uninstall', '--mode', 'unattended']
    cmd_clean_bash = '''grep -v 'export PATH="$PATH:/programs/pathway-tools:"' ~/.bashrc > ~/temp.bashrc; mv ~/temp.bashrc ~/.bashrc'''

    print(' '.join(cmd_uninstall))
    subprocess.call(cmd_uninstall)

    print(cmd_clean_bash)
    subprocess.call(cmd_clean_bash, shell=True)

    if os.path.isdir('/root/AIC-prefs'):
        shutil.rmtree('/root/AIC-prefs')

    ask_delete_ptools(ptools_path)
    return


if __name__ == "__main__":
    main()
