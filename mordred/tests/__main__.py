import os
import nose

if __name__ == '__main__':
    os.chdir(os.path.join(
        os.path.dirname(__file__),
        '..', '..'
    ))
    nose.main(
        defaultTest=os.path.dirname(__file__)
    )
