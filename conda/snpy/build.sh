# the main event
if [ "$(uname)" = "Darwin" ] ; then
   LDFLAGS="$LDFLAGS -undefined dynamic_lookup -bundle"
else
   LDFLAGS="-nostartfiles -shared"
fi
$PYTHON setup.py install
