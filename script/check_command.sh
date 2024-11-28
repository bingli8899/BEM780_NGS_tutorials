#!/bin/bash
# This script checks the dependencies of the software to be used. 

check_dependency() {
    local cmd=$1

    if command -v "$cmd" &>/dev/null; then
        echo -n "$cmd is installed. Version: "
        # Try to get the version of the command, if available
        if $cmd --version &>/dev/null; then
            $cmd --version 2>&1 | head -n 1
        elif $cmd -v &>/dev/null; then
            $cmd -v 2>&1 | head -n 1
        else
            echo "Version information not available."
        fi
    else
        echo "$cmd is not installed."
    fi
}

cmd=$1

check_dependency "$cmd"

