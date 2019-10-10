python2.7 setup.py install
python setup.py install
python2.7 setup.py sdist
python setup.py sdist
twine upload --skip-existing dist/*
