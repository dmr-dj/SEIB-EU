#!/usr/bin/env bash

# Copyright [2023] [Didier M. Roche]
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

set -Eeuo pipefail
trap cleanup SIGINT SIGTERM ERR EXIT

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)

usage() {
  cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") [-h] [-v] -l name_run arg1 [arg2...]

Script description here.

Available options:

-h, --help      Print this help and exit
-v, --verbose   Print script debug info
-l, --label     Name of the run to build
EOF
  exit
}

cleanup() {
  trap - SIGINT SIGTERM ERR EXIT
  # script cleanup here
}

setup_colors() {
  if [[ -t 2 ]] && [[ -z "${NO_COLOR-}" ]] && [[ "${TERM-}" != "dumb" ]]; then
    NOFORMAT='\033[0m' RED='\033[0;31m' GREEN='\033[0;32m' ORANGE='\033[0;33m' BLUE='\033[0;34m' PURPLE='\033[0;35m' CYAN='\033[0;36m' YELLOW='\033[1;33m'
  else
    NOFORMAT='' RED='' GREEN='' ORANGE='' BLUE='' PURPLE='' CYAN='' YELLOW=''
  fi
}

msg() {
  echo >&2 -e "${1-}"
}

die() {
  local msg=$1
  local code=${2-1} # default exit status 1
  msg "$msg"
  exit "$code"
}

parse_params() {
	
  # default values of variables set from params
  outpdataDIR=''
  basedataDIR=''
  nb_years=50
  namerun='DEFAULTSIMU' 
 
  args=("$@")
  
 # [[ ${#args[@]} -lt 1 ]] && die "Missing script arguments: at the minimum I need the label of the run to setup"

  while :; do
    case "${1-}" in
    -h | --help) usage ;;
    -v | --verbose) set -x ;;
    --no-color) NO_COLOR=1 ;;
#    -o | --outputdir) 
#      outpdataDIR="${2-}"
#      shift
#      ;; 
    -n | --nbyears)
      nb_years="${2-}"
      shift
      ;;   
    -l | --label)
      namerun="${2-}"
      shift
      ;;
    -?*) die "Unknown option: $1" ;;
    *) break ;;
    esac
    shift
  done

  # check required arguments
#  [[ -z "${basedataDIR-}" ]] && die "Missing required directory: inputdata DIR"
#  [[ -z "${outpdataDIR-}" ]] && die "Missing required directory: outptdata DIR"
#  [[ ! -d "${basedataDIR-}" ]] && die "Inputdata argument must be an existing directory"
#  [[ ! -d "${outpdataDIR-}" ]] && die "Outptdata argument must be an existing directory"

  return 0
}

parse_params "$@"

setup_colors

mkdir -p wkdir/${namerun}
cp -p src/* wkdir/${namerun}/.
cp -p params/* wkdir/${namerun}/.
sed -e "s%THE_NUMBER_OF_YEARS_TO_RUN%${nb_years}%g" params/parameter.txt > wkdir/${namerun}/parameter.txt
cp -p config/Makefile wkdir/${namerun}/.
cp -pr visualizer/ wkdir/${namerun}/.
mkdir wkdir/${namerun}/result_tmp
mkdir wkdir/${namerun}/result_visualize
cp -p SEIB-IAB.job wkdir/${namerun}/.

msg "${GREEN}Finalized setup in wkdir/${namerun}${NOFORMAT}"

# The End of All Things (op. cit.)

