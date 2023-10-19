# Some tests

print("Testing GS unit computation")
print("-" * 20)
print("\nD = 24, p = 7")

if str(GS_unit(24, 7)) == "7*x^2 + 13*x + 7":
    print("PASSED")
else:
    print("FAILED")

print("\nD = 205, p = 2")
if str(GS_unit(205, 2, 30, 40)) == "4*x^4 - 5*x^3 + 7*x^2 - 5*x + 4":
    print("PASSED")
else:
    print("FAILED")

print("\nD = 321, p = 7")
if str(
        GS_unit(321, 7, 40, 50)
) == "4747561509943*x^6 - 9414439453335*x^5 + 6095965902306*x^4 - 1158716919139*x^3 + 6095965902306*x^2 - 9414439453335*x + 4747561509943":
    print("PASSED")
else:
    print("FAILED")

print("-"*20 + "\n")
print("Testing diagonal restriction computation")
print("-"*20 + "\n")

print("D=69, k in {1,..,4}")
if str([diagonal_restriction(BQFs(69)[0],k) for k in [1..4]]) ==  '[0, 2 + 480*q + 4320*q^2 + 13440*q^3 + 35040*q^4 + 60480*q^5 + O(q^6), -32/3 + 5376*q + 177408*q^2 + 1311744*q^3 + 5682432*q^4 + 16805376*q^5 + O(q^6), 165 + 79200*q + 10216800*q^2 + 173289600*q^3 + 1307829600*q^4 + 6187579200*q^5 + O(q^6)]':
    print("PASSED")
else:
    print("FAILED")

print("\nD=24*7^2, k in {1,..,4}")
if str([diagonal_restriction(BQFs(24*7^2)[0],k) for k in [1..4]]) ==  '[-2 - 8*q - 24*q^2 - 32*q^3 - 56*q^4 - 48*q^5 + O(q^6), 9668*q + 85052*q^2 + 264824*q^3 + 690084*q^4 + 1196608*q^5 + O(q^6), 479780*q + 15506204*q^2 + 114357992*q^3 + 495641076*q^4 + 1465429168*q^5 + O(q^6), 102280844*q + 13240347284*q^2 + 224482262792*q^3 + 1694155192044*q^4 + 8015255739616*q^5 + O(q^6)]':
    print("PASSED")
else:
    print("FAILED")

print("\nD=105, k = 1, 11-stabilised")

if diagonal_restriction(BQFs(105)[0],1,11).q_expansion() == 0:
    print("PASSED")
else:
    print("FAILED")
