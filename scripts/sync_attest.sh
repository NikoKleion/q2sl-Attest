#!/bin/sh
# copy the committed attest/ subfolder into the qrng-attest repository and push it
set -e
repo=https://github.com/NikoKleion/qrng-attest.git
rev=$(git rev-parse --short HEAD)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
git clone --quiet "$repo" "$work"
git -C "$work" rm -rq --ignore-unmatch .
git archive HEAD:attest | tar -x -C "$work"
git -C "$work" add -A
if git -C "$work" diff --cached --quiet; then
    echo "qrng-attest already matches attest/ at $rev"
    exit 0
fi
git -C "$work" commit --quiet -m "attest/ from q2sl attest $rev"
git -C "$work" push --quiet origin HEAD
echo "pushed attest/ at $rev to $repo"
