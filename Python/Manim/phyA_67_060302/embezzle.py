from manim import *
import numpy as np

config.frame_width = 12
config.pixel_width = 1280
config.pixel_height = 1024

class ExpandingEmbezzlingState(Scene):
    def construct(self):
        # 10개 막대 초기 상태 (이전 장면의 마지막 상태와 동일)
        initial_10_values = [1/np.sqrt(i) for i in range(1, 11)]
        weight_10 = np.sum([1/i for i in range(1, 11)])
        initial_10_values = [i/np.sqrt(weight_10) for i in initial_10_values]
        
        # 20개 막대 확장 상태 - np.interp 사용
        x_old = np.arange(1, 11)  # 기존 10개 위치 (1부터 10까지)
        x_new = np.linspace(1, 10, 20)  # 새로운 20개 위치 (1부터 10까지를 20등분)

        # 기존 10개 값에 대한 보간으로 20개 값 생성
        expanded_20_values = np.interp(x_new, x_old, initial_10_values).tolist()
        
        max_height = 2
        bar_width_10 = 0.2
        bar_width_20 = 0.15  # 20개로 늘어나면서 막대 폭 줄임
        spacing_10 = 0.3
        spacing_20 = 0.2  # 간격도 줄임
        
        # 10개 막대 생성 (빨간색으로 시작)
        bars_10 = VGroup()
        labels_10 = VGroup()
        
        for i in range(10):
            height = initial_10_values[i] * max_height / 0.5
            bar = Rectangle(
                width=bar_width_10,
                height=height,
                fill_color=RED,
                fill_opacity=0.7,
                stroke_width=2
            )
            bar.shift(RIGHT * (i - 4.5) * spacing_10)
            bar.shift(UP * height/2)
            
            label = Text(f"j_{i+1}", font_size=12)
            label.next_to(bar, DOWN, buff=0.1)
            
            bars_10.add(bar)
            labels_10.add(label)
        
        # 20개 막대 생성 (목표 상태)
        bars_20 = VGroup()
        labels_20 = VGroup()
        
        for i in range(20):
            height = expanded_20_values[i] * max_height / 0.5
            bar = Rectangle(
                width=bar_width_20,
                height=height,
                fill_color=TEAL,
                fill_opacity=0.7,
                stroke_width=2
            )
            bar.shift(RIGHT * (i - 9.5) * spacing_20)
            bar.shift(UP * height/2)
            
            label = Text(f"j_{i+1}", font_size=10)
            label.next_to(bar, DOWN, buff=0.1)
            
            bars_20.add(bar)
            labels_20.add(label)
        
        # 제목들
        title1 = Text("Embezzling State (n=10)", font_size=18)
        title1.to_edge(UP)
        
        title2 = Text("Expanded Embezzling State (n=20)", font_size=18)
        title2.to_edge(UP)
        
        # 애니메이션 시퀀스
        
        # 부드러운 전환을 위해 중간 단계 추가
        self.play(Write(title1))
        self.play(Create(bars_10))
        self.play(Write(labels_10))
        self.wait(2)
        
        # 막대들을 중앙으로 압축
        self.play(
            bars_10.animate.scale(0.7),
            labels_10.animate.scale(0.7)
        )
        
        # 제목 변경
        self.play(Transform(title1, title2))
        
        # 추가 막대들을 양쪽에서 나타나게 하기
        left_bars = VGroup(*bars_20[:5])   # 왼쪽 5개
        middle_bars = VGroup(*bars_20[5:15])  # 중간 10개 (기존 위치)
        right_bars = VGroup(*bars_20[15:])    # 오른쪽 5개
        
        # 기존 막대를 중간으로 변환
        self.play(Transform(bars_10, middle_bars))
        
        # 양쪽에 새 막대 추가
        self.play(
            FadeIn(left_bars, shift=RIGHT),
            FadeIn(right_bars, shift=LEFT)
        )
        
        # 라벨 업데이트
        self.play(Transform(labels_10, labels_20))
        
        self.wait(3)
