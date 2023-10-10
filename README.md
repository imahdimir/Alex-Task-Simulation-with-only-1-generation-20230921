# Apple ARM processor machines

There can be difficulties install *snipar* on Apple ARM processor machines due
to lack of available versions of scientific computing software made for these
processors' architectures. A workaround for this is to use *Snipar* in a docker
container.

The following steps will guide you on how to create a suitable container for
your M1/M2 MacBook and seamlessly integrate it into a VSCode environment. These
steps assume you have little knowledge about Docker, so we will start from
scratch.

## 1. Installing Docker Engine

Ensure that Docker Engine is installed on your machine. If it's already
installed, you can skip this step. Otherwise, you can install it by following
the instructions provided at the following link:

[Install Docker Engine on a MacBook](https://docs.docker.com/desktop/install/mac-install/)

## 2. Creating the Docker Container

To install the Snipar package, you need to create a Docker container that
emulates a virtual Ubuntu machine where Snipar can be installed.

To create the appropriate Docker container, follow these steps:

- Start the Docker engine. Just open the Docker application from your
  Applications folder.

- While the engine is running, open your terminal and execute the following
  command to create a container named "snipar_container" (you can choose a
  different name if you prefer):

  ```bash
  docker run --name snipar_container -it amd64/python:3.9.9-slim-buster /bin/bash
  ```

After running this command, you should see the "snipar_container" listed in the
Docker GUI under the "Containers" tab.

- You can close the terminal at this point.

## 3. Running the Container

To use the environment of the created container within VSCode or other IDEs,
ensure that the container is running. You can do this easily through the Docker
GUI:

- Open the Docker GUI.
- In the "Containers" section of the dashboard, locate "snipar_container" (or
  the name you chose for your container).
- Click the play button to start the container.

Once the container is running, you can access its terminal and files through the
Docker GUI. Keep in mind that the files within the container are isolated from
your macOS files. You won't be able to access them using Finder, but you can
manage them through the Docker GUI, including uploading and downloading.

## 4. Attaching the Container to VSCode

If you prefer a smoother development experience with VSCode, follow these steps
to attach the container to VSCode:

- Install the "Dev Containers" extension by visiting the following link in the
  VS Marketplace:

  [Dev Containers Extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

- Once installed, a new icon will appear on the left sidebar of VSCode labeled "
  Remote Explorer."

- Click on "Remote Explorer" and, under the "Dev Containers" section, locate "
  snipar_container."

- Next to the container name, you'll find a button that says "Attach in a new
  window." Click on this button, and VSCode will open a new window attached to
  the container.

In the newly attached VSCode window, the terminal will be connected to the
container, similar to the "Exec" tab in the Docker GUI. You can also open
folders from the container environment (not your Mac itself) using the "Open
Folder" button in VSCode. This makes it more convenient to manage files while
writing code, and you can run your modules using the terminal directly within
VSCode.

## 5. Installing *SNIPAR* Package in the Container

Up to this point, you've set up the environment for Snipar but haven't installed
it in the container. To install Snipar, follow these steps:

- Open the terminal in the VSCode window attached to "snipar_container" (or use
  the "Exec" tab in the Docker GUI).
- Run the following commands:
  ```bash
  pip install --upgrade pip
  pip install snipar
  ```

If you're conducting the simulation exercise of *Snipar* it requires the "plink"
package, we recommend installing it natively on your Mac and using the
simulation results. You can download the simulation result files from the
container to your macOS environment and perform the analysis without the need
for a Docker container.
