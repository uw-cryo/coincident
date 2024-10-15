#!/bin/bash

# Check the shell being used
if [[ "$SHELL" == *"/bash" ]]; then
    echo "Running activation for Git Bash..."
    export PATH="$PATH:/c/Program Files/Git/bin"

elif [[ "$SHELL" == *"/zsh" ]]; then
    echo "Running activation for Zsh..."
    export PATH="$PATH:/usr/bin/zsh"  # Adjust as necessary

else
    echo "Running activation for other Unix-like shells..."
    export PATH="$PATH:/usr/local/bin"  # Adjust as necessary
fi

# Common commands for all shells
echo "Activating the project environment..."
# Any additional setup can go here, e.g., setting virtualenv paths, etc.