#!/bin/bash

#Парметры ниже можно определить явно привызове скрипта либо через # в самом скрипте
#Определяет рабочую директорию для задания. Если не задана, то рабочей является текущая директория:

# Эта папка должна существовать 
#SBATCH -D /mnt/pool/rhic/1/fkurushin/McPicoReader/Resolution/TMP

#Определяет имя задачи:
#SBATCH -J readTree

#Задаёт очередь в которую добавляется задача, если SBATCH -p queue, что делает строчка снизу - неизвестно: 
#SBATCH -p compute

#Запрашивает N вычислительных узлов для задачи:
#SBATCH -N 1

#Запрашивает N процессов для задачи:
#SBATCH -n 1

#Ограничивает максимальное время выполнения задачи:
#SBATCH --time=0-00:20:00

#SBATCH --array=1-39

#Эти папки тоже должны существовать
#Задаются имена файлов ошибок (stderr), и стандартного вывода (stdout).
#SBATCH -o /mnt/pool/rhic/1/fkurushin/McPicoReader/Resolution/SGE_OUT/slurm_%A_%a.out
#SBATCH -e /mnt/pool/rhic/1/fkurushin/McPicoReader/Resolution/SGE_OUT/slurm_%A_%a.err


#This will be set to the first job ID of the array.
export JOB_ID=$SLURM_ARRAY_JOB_ID
#This will be set to the job array index value.
export TASK_ID=$SLURM_ARRAY_TASK_ID


export MAIN_DIR=/mnt/pool/rhic/1/fkurushin/McPicoReader/Resolution
export FILELIST_DIR=${MAIN_DIR}/lists

#Команды ниже копируют аутпут после выполнения команд в кавычках 
IN_LIST=`ls $FILELIST_DIR`
#Получаем лист с номером task id
INPUT=`ls $FILELIST_DIR | sed "${TASK_ID}q;d"`

export OUT=$MAIN_DIR/OUT
export OUT_LOG=$OUT/log
export OUT_FILE=$OUT/files

export LOG=$OUT_LOG/JOB_${JOB_ID}_${TASK_ID}.log
export OUTPUT=$OUT_FILE/JOB_${JOB_ID}_${TASK_ID}.root

mkdir -p $OUT_LOG
mkdir -p $OUT_FILE

export ROOT_MACRO=${MAIN_DIR}/readTree.C

source /mnt/pool/nica/5/parfenovpeter/Soft/Cherenkov/ROOT/build-cxx11/bin/thisroot.sh

root -l -b -q $ROOT_MACRO'("'$INPUT'","'$OUTPUT'")' &>> $LOG

echo "Job finished!" &>> $LOG
