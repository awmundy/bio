import datetime
import os
import sys
sys.path.append('/home/awmundy/code/')
from airflow import DAG
from airflow.operators.python import PythonOperator
from bio.data_prep.config import cfg, run

default_args = {
    'owner': 'awmundy',
    'retries': 5,
    'retry_delay': datetime.timedelta(minutes=1),
    }

cfg = {'ref_genome': ''}
threads = 2

def kallisto_build_index(ref_genome):
    print(f'index constructed with {ref_genome}')

def prep_fastq_files(cfg, run):
    print('fastq files prepped')

def fastqc(rna_txs_dir, threads):
    print('fastqc done')

def kallisto_quant(index_fpath, rna_txs_dir, threads, seq_params):
    print('kallisto quant done')

def multiqc(rna_txs_dir):
    print('multiqc done')


# index_fpath = kallisto_build_index(cfg['ref_genome'])
# prep_fastq_files(cfg, run)
# fastqc(cfg['rna_txs_dir'], threads)
# kallisto_quant(index_fpath, cfg['rna_txs_dir'], threads, cfg['seq_params'])
# multiqc(cfg['rna_txs_dir'])

with DAG(
        default_args=default_args,
        dag_id='airflow_dag',
        description='testing dags',
        start_date=datetime.datetime(2022, 2, 2),
        schedule_interval='@once',
        ) as dag:
    task1 = PythonOperator(
            task_id='kallisto_build_index',
            python_callable=kallisto_build_index,
            op_kwargs={'ref_genome': cfg['ref_genome']},
            )
    # task2 = PythonOperator(
    #         task_id='downstream_task_1',
    #         python_callable=task_2,
    #         )
    # task3 = PythonOperator(
    #         task_id='downstream_task_2',
    #         python_callable=task_3,
    #         )
    #
    # task1.set_downstream(task2)
    # task1.set_downstream(task3)
    #
    # task1 >> [task2, task3]
    task1