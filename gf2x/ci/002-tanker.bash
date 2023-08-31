
exported_variables=(
      CI_JOB_NAME
      CI_COMMIT_SHORT_SHA
      CI_JOB_ID
      CI_JOB_STAGE
      CI_PROJECT_NAMESPACE
      CI_PROJECT_NAME
      DOCKER_SCRIPT
)

exports=(RUNTIME_TYPE=tanker)
for v in "${exported_variables[@]}" ; do
    if [ "${!v}" ] ; then
        exports+=("$v=\"${!v}\"")
    fi
done

tanker() {
    "$(dirname $0)/utilities/tanker/tanker.sh" -B "${TANKER_DATABASE:-/tmp}/tanker-$(id -n -u)" "$@"
}

