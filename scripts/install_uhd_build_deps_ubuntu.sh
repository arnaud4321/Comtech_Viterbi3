#!/usr/bin/env bash
# Ubuntu: packages recommended to *build* UHD from source (Ettus manual, "Setting up the dependencies on Ubuntu").
# https://files.ettus.com/manual/page_build_guide.html
set -euo pipefail
exec sudo apt-get install -y \
  autoconf automake build-essential ccache cmake cpufrequtils doxygen ethtool \
  g++ git inetutils-tools libboost-all-dev libncurses5 libncurses5-dev \
  libusb-1.0-0 libusb-1.0-0-dev libusb-dev python3-dev python3-mako python3-numpy \
  python3-requests python3-scipy python3-setuptools python3-ruamel.yaml
