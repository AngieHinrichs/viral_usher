# 'build' subcommand: build a tree according to parameters in config file

import datetime
import filecmp
import importlib.metadata
import os
import shutil
import subprocess
import sys
from . import config

docker_platform = 'linux/amd64'
docker_workdir = '/data'


def choose_container_runtime():
    """Detect available container runtime. Prefer Apptainer if available.
    Returns a tuple (runtime, error_message). runtime is one of 'apptainer', 'docker' or None on failure."""
    # Check for apptainer first (preferred)
    try:
        subprocess.run(['apptainer', '--version'], check=True, capture_output=True)
        return 'apptainer', None
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    # Fallback to docker CLI
    try:
        subprocess.run(['docker', '--help'], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None, "Unable to run 'docker' or 'apptainer' -- please install one of them and try again."
    try:
        subprocess.run(['docker', 'images'], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        return None, "Unable to run 'docker images' -- please make sure the docker daemon is running (e.g. on a Mac, start the Docker app)"
    return 'docker', None


def pull_docker_image_subprocess(image_name: str):
    """Pull a docker image using the docker CLI."""
    try:
        subprocess.run(['docker', 'pull', '--platform', docker_platform, image_name], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to pull docker image {image_name}:\n{e}", file=sys.stderr)
        sys.exit(1)
    except TimeoutError:
        print("Timeout error while trying to pull docker image; maybe try again later?")
        sys.exit(1)


def get_container_image(runtime: str, args_docker_image: str):
    """Return an image reference appropriate for the runtime. For docker runtime, ensure the image exists locally (pull if needed).
    For apptainer runtime, return a docker://... URI so apptainer can pull it on demand."""
    if args_docker_image:
        image_name = args_docker_image
    else:
        viral_usher_version = importlib.metadata.version('viral_usher')
        image_name = config.DEFAULT_DOCKER_IMAGE + ':v' + viral_usher_version

    if runtime == 'docker':
        # Check if image exists; if not, pull it
        try:
            subprocess.run(['docker', 'image', 'inspect', image_name], check=True, capture_output=True)
            print(f"Using docker image {image_name}")
        except subprocess.CalledProcessError:
            print(f"Docker image {image_name} is not present, pulling it from DockerHub...")
            pull_docker_image_subprocess(image_name)
        return image_name
    else:  # apptainer
        # apptainer accepts docker://<image> URIs; ensure we return that form
        if '://' in image_name:
            return image_name
        return f"docker://{image_name}"


def maybe_copy_to_workdir(filepath, workdir):
    """If filepath is a file not an URL, return filepath's relative path within workdir, copying the file into workdir if necessary."""
    if filepath.startswith('http://') or filepath.startswith('https://'):
        return filepath
    basename = os.path.basename(filepath)
    dirname = os.path.dirname(filepath)
    dirname_abs = os.path.abspath(dirname)
    workdir_abs = os.path.abspath(workdir)
    if dirname_abs == workdir_abs:
        return basename
    else:
        if dirname_abs.startswith(workdir_abs + '/'):
            relpath = dirname_abs.removeprefix(workdir_abs + '/')
            relname = relpath + '/' + basename
            return relname
        else:
            workdir_file = workdir + '/' + basename
            if os.path.exists(workdir_file):
                if not filecmp.cmp(filepath, workdir_file, shallow=False):
                    print(f"File {filepath} is specified, but workdir file {workdir_file} already exists and is not the same", file=sys.stderr)
                    sys.exit(1)
            else:
                shutil.copy(filepath, workdir)
            return basename


def rewrite_config(config_contents, workdir, no_genbank):
    """Make a new timestamped config file in workdir that contains the given config settings"""
    now = datetime.datetime.now()
    new_name = 'local_config.' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.toml'
    new_path = workdir + '/' + new_name
    config.write_config(config_contents, new_path, check_paths=False, no_genbank=no_genbank)
    return new_name


def localize_file(config_contents, key, workdir):
    """If key is in config_contents, check whether the file is in workdir; if not, copy it there and update config_contents[key].
    Return True if config_contents was changed, else False."""
    if key in config_contents and config_contents[key] != "":
        filepath = config_contents[key]
        filepath_rel = maybe_copy_to_workdir(filepath, workdir)
        if filepath_rel != filepath:
            config_contents[key] = filepath_rel
            return True
    return False


def localize_config(args_config, config_contents, no_genbank):
    """When running in docker, only the workdir will be available (as the current directory).
    So absolute paths in the config need to be converted to relative paths, and if input files
    are not already in the workdir, then they need to be copied there."""
    workdir = config_contents['workdir']
    # If any user-provided files are given, make sure they will be available in workdir and that
    # the config used by the docker image script uses their relative paths in workdir.
    config_changed = False
    for local_file_setting in ['extra_fasta', 'extra_metadata', 'ref_fasta', 'ref_gbff', 'taxonium_overlay_html']:
        config_changed |= localize_file(config_contents, local_file_setting, workdir)
    if config_changed:
        config_rel = rewrite_config(config_contents, workdir, no_genbank)
    else:
        config_rel = maybe_copy_to_workdir(args_config, workdir)
    return config_rel


def run_in_container(runtime, workdir, config_rel, args_docker_image, args_update, args_no_genbank, args_threads):
    """Run the pipeline using either docker CLI or apptainer CLI depending on runtime."""
    user_id = os.getuid()
    group_id = os.getgid()
    image_ref = get_container_image(runtime, args_docker_image)
    base_command = ["viral_usher_build", "--config", config_rel]
    if args_update:
        base_command.append("--update")
    if args_no_genbank:
        base_command.append("--no_genbank")
    if args_threads:
        base_command += ["--threads", str(args_threads)]

    if runtime == 'docker':
        cmd = [
            'docker', 'run', '--rm', '--platform', docker_platform,
            '-w', docker_workdir,
            '--network', 'host', '--user', f"{user_id}:{group_id}",
            '-v', f"{workdir}:{docker_workdir}",
        ] + [image_ref] + base_command
    else:  # apptainer
        # apptainer will pull docker:// images automatically; bind workdir
        # Use 'exec' to run the command inside the image
        cmd = [
            'apptainer', 'exec', '--bind', f"{workdir}:{docker_workdir}",
            '--pwd', docker_workdir,
            image_ref
        ] + base_command

    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            print(line, end='')
        proc.wait()
        if proc.returncode != 0:
            print(f"{runtime} container {image_ref} failed with exit code {proc.returncode}", file=sys.stderr)
            sys.exit(proc.returncode)
    except FileNotFoundError as e:
        print(f"Failed to run container runtime: {e}", file=sys.stderr)
        sys.exit(1)
    except TimeoutError:
        print(f"Timeout error while trying to run {runtime}; maybe try again later?")
        sys.exit(1)


def handle_build(args):
    """Set up input files in workdir, run pipeline in docker and check final results."""
    runtime, error = choose_container_runtime()
    if runtime is None:
        print(f"\n{error}\n", file=sys.stderr)
        sys.exit(1)
    print(f"Running pipeline with config file: {args.config}")
    try:
        # Don't download URLs to local temp files here; let the docker container download & use them
        config_contents = config.parse_config(args.config, resolve_url_keys=False, no_genbank=args.no_genbank)
    except (FileNotFoundError, PermissionError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    config_rel = localize_config(args.config, config_contents, args.no_genbank)
    workdir = config_contents['workdir']
    run_in_container(runtime, workdir, config_rel, args.docker_image, args.update, args.no_genbank, args.threads)
    if os.path.exists(f"{workdir}/tree.jsonl.gz"):
        print(f"Success -- you can view {workdir}/tree.jsonl.gz using https://taxonium.org/ now.")
    else:
        print(f"The pipeline exited without an error status but the expected result file {workdir}/tree.jsonl.gz is not there.  Please file an issue in GitHub: https://github.com/AngieHinrichs/viral_usher/issues/new", file=sys.stderr)
        sys.exit(1)
