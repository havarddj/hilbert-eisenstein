import csv


def batch_compute_D(upper_bound, p, lower_bound=1, replace=False):

    # first file create if it doesn't exist
    with open(f"data/D{upper_bound}p{p}.csv", 'w+', newline='') as csvfile:
        keys = ["D", "P", "p"]
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        if not writer:
            writer.writeheader()

    D = lower_bound
    while True:
        # check if D is a good discriminant:
        while not is_GS_disc(D, p):
            D += 1

        if D > upper_bound:
            return 0

        # check if already computed:
        if not replace:
            with open(f"data/D{upper_bound}p{p}.csv", 'r',
                      newline='') as csvfile:
                keys = ["D", "P", "p"]
                reader = csv.DictReader(csvfile, fieldnames=keys)

                if D in reader:
                    if D["P"] != 0:
                        D += 1

        P = 0
        pprec = 30
        nterms = min(
            ModularForms(weight=2 +
                         (p - 1) * floor(pprec * (p + 1) / p)).dimension() - 1,
            pprec)
        while P == 0:
            P = GS_unit(D, p, pprec=pprec, nterms=pprec + 20)
            pprec += 30
            nterms += 30

        with open(f"data/D{upper_bound}p{p}.csv", 'a+', newline='') as csvfile:
            keys = ["D", "P", "p"]
            writer = csv.DictWriter(csvfile, fieldnames=keys)

            writer.writerow({'D': D, 'P': P, 'p': p})
        D += 1
