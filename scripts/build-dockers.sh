#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
pushd "$repo_dir" &> /dev/null

tag=""

function usage {
    cat << EOF
Usage: $0 [ARGS] -- [projects]

Build project dockers. If no project is specified, build all.
Update .env with the list of current dockers.

ARGS:
    -h|--help: Show this message and exit.
    -t|--tag: Specify docker tag. Defaults to the latest reachable Git tag.
    --push/--no-push: Whether or not to push images to the remote repo
EOF
}

args=()
projects=()
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h"|"--help")
            usage
            exit 0
            ;;
        "-t"|"--tag")
            if [[ "$#" -lt 2 || -z "$2" || "$2" == -* ]]; then
                echo "--tag requires a tag value." >&2
                exit 2
            fi
            tag="$2"
            shift 2
            ;;
        "--"?*)
            args+=("$1")
            shift 1
            ;;
        "--")
            shift 1
            break
            ;;
        *)
            break
            ;;
    esac
done

if [[ "$#" -gt 0 ]]; then
    projects+=("${@}")
fi

if [[ -z "$tag" ]]; then
    if ! tag=$(git describe --tags --abbrev=0 2>/dev/null); then
        echo "No Git tag found; pass --tag TAG." >&2
        exit 2
    fi
fi

case ${#projects[@]} in
    0)
        all_projects="true"
        ;;
    1)
        if [[ "${projects[0]}" == "all" ]]; then
            all_projects="true"
        else
            all_projects="false"
        fi
        ;;
    *)
        all_projects="false"
        ;;
esac

tmp_env=$(mktemp)
trap 'rm -rf "$tmp_env"' EXIT

"$script_dir/find-projects.sh" \
| while read -r project type; do
    if [[ "$all_projects" == "false" ]] && [[ ! " ${projects[*]} " =~ [[:space:]]${project}[[:space:]] ]]; then
        continue
    fi
    case "$type" in
        "pixi")
            1>&2 echo "Building docker for pixi project: $project:$tag"
            "$script_dir/build-project-docker.sh" ${args[@]+"${args[@]}"} "$project" "$tag"
            ;;
        "uv")
            1>&2 echo "Building docker for uv project: $project:$tag"
            "$script_dir/build-project-docker.sh" ${args[@]+"${args[@]}"} "$project" "$tag"
            ;;
        "yaml")
            1>&2 echo "Building docker for yaml project: $project:$tag"
            "$script_dir/build-yaml-docker.sh" ${args[@]+"${args[@]}"} "$project" "$tag"
            ;;
        *)
            1>&2 echo "Unknown project type: $type for project: $project"
            exit 1
            ;;
    esac
done \
| awk -F = \
    '
    { envs_dict[$1]=$2 }
    END {
        for(env in envs_dict) {
            printf "%s=%s\n", env, envs_dict[env]
        }
    }
    ' \
    "$repo_dir/.env" \
    - \
> "$tmp_env"

mv "$tmp_env" "$repo_dir/.env"
