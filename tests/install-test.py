"""
This file does a very quick check to see if you have installed everything
"""

def main():
    """
    Try to import
    """
    try:
        from humbert import arithmetic
        from sympy import randprime
        p = randprime(10 ** 3, 10 ** 4)
        symb = arithmetic.kronecker(3, p)
        print(f"The Legendre symbol (3|{p}) is {symb}")
        print("Everything went smoothly it seems")
    except ImportError:
        print("We couldn't import the humbert package...")
        print("Did you run `python3 -m pip install .` ?")

    except Exception as e:
        print(f"Unexpected error! {e}")
        print("Disaster, please contact me [sam.frengley@bristol.ac.uk]")


if __name__ == "__main__":
    main()
