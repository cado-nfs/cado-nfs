
exported_variables=(
      BUILD_NAME
      COMMIT_SHORT_SHA
      JOB_ID
      REPOSITORY
      DOCKER_SCRIPT
)

exports=(RUNTIME_TYPE=tanker)
for v in "${exported_variables[@]}" ; do
    if [ "${!v}" ] ; then
        exports+=("$v=\"${!v}\"")
    fi
done

tanker() {
    # add this because we want to always have a current list of images
    (cd "$(dirname $0)/utilities/tanker" ; rm -f images.txt || : ; make images.txt)
    "$(dirname $0)/utilities/tanker/tanker.sh" -B "${TANKER_DATABASE:-/tmp}/tanker-$(id -n -u)" "$@"
}

