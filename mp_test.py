from multiprocessing import Pool


def main():
    for i in range(0,100000):
        a = i**i

if __name__ == '__main__':
    main()