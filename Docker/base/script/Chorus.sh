#!/bin/bash

echo $CHORUS_USER $CHORUS_HOME

BLUE='\033[0;34m'
GREEN='\033[0;32m'
NC='\033[0m'

if ! id -u $CHORUS_USER > /dev/null 2>&1; then
  echo -e "use local user:  ${BLUE}$CHORUS_USER${NC}"
  addgroup --gid $CHORUS_UID  $CHORUS_USER
  adduser --force-badname --disabled-password  --gecos "" --uid  $CHORUS_UID --gid $CHORUS_UID $CHORUS_USER
fi

if [ -e $CHORUS_HOME ];
    then
      echo -e "${GREEN}$CHORUS_HOME${NC} exists"
      chown -R $CHORUS_USER $CHORUS_HOME
else
    mkdir -p $CHORUS_HOME
    chown -R $CHORUS_UID $CHORUS_HOME
fi

for arg in "$*"
do
  	exec python3 /opt/software/Chorus/Chorus.py --docker True $*
done
