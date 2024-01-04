ALL_SPACES = list('123456789')
X,O,BLANK = 'X','O',' '

def main():
    print('틱택토  게임에 오신 당신을 환영합니다!')
    gameBoard = getBlankBoard()
    currentPlayer, nextPlayer = X, O

    while true:
        print(getBoardsStr(gameBoard))

        move = 0
        while not isValidSpace(gameboard,move):
            print(f'{currentPlayer}의 움직임은?(1-9)')
            move = input()
        updateBoard(gameBoard, move, currentPlayer)

        if isWInner(gameBoard, currentPlayer):
            print(gettBoardStr(gameBoard))
            print(currenyPlayer + '가 승리했습니다!')
            break
        elif isBoardfull(gameBoard):
            print(getBoardStr(gameBoard))
            print('무승부 게임입니다!')
            break
        currentPlayer, nextPlayer
    print('즐겁게 퍼즐을 풀어주셔서 감사합니다!')

def getBlankBoard():
        """비어있는 새 틱택토 말판을 생성한다."""
        return f'''
        {board['1']}|{board['2']}|{board['3']} 1 2 3
        -+-+-
        {board['4']}|{board['5']}|{board['6']} 4 5 6
        -+-+-
        {board['7']}|{board['8']}|{board['9']} 7 8 9 '''

def isValidSpace(board, space):
     """boarddml .space가 유효한 칸 번호이며, 그 칸이 비어 있을 경우 True를 반환한다."""
     return 0 < int(space) < 10 and (space in ALL_SPACES or board[space] == BLANK)

def isWinner(board, player):
    """player가 이 board에서 승자인 경우 True를 반환한다."""
    b, p = board, player
    return ((b['1'] == b['2'] == b['3'] == p) or
            (b['4'] == b['5'] == b['6'] == p) or
            (b['7'] == b['8'] == b['9'] == p) or
            (b['1'] == b['4'] == b['7'] == p) or
            (b['2'] == b['5'] == b['8'] == p) or
            (b['3'] == b['6'] == b['9'] == p) or
            (b['3'] == b['5'] == b['7'] == p)or
            (b['1'] == b['5'] == b['9'] == p))

def isBoardFull(board):
    """board의 모든 칸이 차 있다면 True를 반환한다."""
    for space in ALL_SPACES:
          if board[space] == BLANK:
               return False
    return True

def updatedBoard(board, space, mark):
    """boarddml space를 mark로 설정한다."""
    board[space] = mark

if __name__=='main':
    main()
