from rs2_score_calculator import get_rs2_score as calculate_rs2

if __name__ == '__main__':
	sequence = 'AAAAAAAAAAATAAAAAAAAAAAAAGGAAA'
	rs2_scores = calculate_rs2(sequence), calculate_rs2(sequence)
	print(rs2_scores)