#!/bin/bash

echo Executing command:
echo $CMD
matlab -nodisplay -nosplash -nojvm -r $"$CMD; exit"
