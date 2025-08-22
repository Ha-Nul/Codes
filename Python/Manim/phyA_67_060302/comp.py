from manim import *
import numpy as np

config.frame_width = 12
config.pixel_width = 1280
config.pixel_height = 1024

class ComparingBars(Scene):
    def construct(self):
        # 10개 막대 초기 상태
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
        bar_width_20 = 0.15
        spacing_10 = 0.3
        spacing_20 = 0.2
        
        # 10개 막대 생성 (빨간색으로 시작)
        bars_10 = VGroup()
        #labels_10 = VGroup()
        
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
            
            #label = Text(f"j_{i+1}", font_size=12)
            #label.next_to(bar, DOWN, buff=0.1)
            
            bars_10.add(bar)
            #labels_10.add(label)
        
        # 20개 막대 생성 (목표 상태)
        bars_20 = VGroup()
        #labels_20 = VGroup()
        
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
            
            #label = Text(f"j{i+1}", font_size=10)
            #label.next_to(bar, DOWN, buff=0.1)
            
            bars_20.add(bar)
            #labels_20.add(label)
        
        # 제목들
        title1 = Text("Embezzling State (n=10)", font_size=18)
        title1.to_edge(UP)
        
        title2 = Text("Comparing 5th Bar with Expanded State", font_size=18)
        title2.to_edge(UP)
        
        title3 = Text("5th Bar Position in Expanded State", font_size=18)
        title3.to_edge(UP)
        
        # 애니메이션 시퀀스
        
        # 1. 초기 10개 막대 표시
        self.play(Write(title1))
        self.play(Create(bars_10))
        #self.play(Write(labels_10))
        self.wait(2)
        
        # 2. 5번째 막대 강조 (인덱스 4)
        fifth_bar = bars_10[4]
        #fifth_label = labels_10[4]
        
        # 5번째 막대를 강조표시
        self.play(
            fifth_bar.animate.set_fill(YELLOW, opacity=1),
            fifth_bar.animate.set_stroke(YELLOW, width=4)
        )
        self.wait(1)
        
        # 3. 제목 변경
        self.play(Transform(title1, title2))
        
        # 4. 5번째 막대와 라벨을 제외한 나머지 요소들 페이드아웃
        other_bars = VGroup(*[bars_10[i] for i in range(10) if i != 4])
        #other_labels = VGroup(*[labels_10[i] for i in range(10) if i != 4])
        
        self.play(
            FadeOut(other_bars),
            #FadeOut(other_labels)
        )
        
        # 5. 20개 막대 페이드인 (5번째 막대는 유지)
        self.play(
            FadeIn(bars_20),
            #FadeIn(labels_20)
        )
        self.wait(1)
        
        # 6. 제목 변경
        self.play(Transform(title1, title3))
        
        # 7. 5번째 막대를 20개 막대의 10번째 위치로 이동 (인덱스 9)
        target_bar = bars_20[4]  # 10번째 막대 (인덱스 9)
        target_position = target_bar.get_center()
        
        # 5번째 막대를 10번째 막대 위치로 이동
        self.play(
            fifth_bar.animate.move_to(target_position),# + LEFT * 0.3),  # 약간 왼쪽으로 오프셋
            #fifth_label.animate.move_to(target_position + DOWN * 0.8) #+ LEFT * 0.3)
        )
        
        # 8. 10번째 막대도 강조하여 비교
        self.play(
            target_bar.animate.set_fill(ORANGE, opacity=1),
            target_bar.animate.set_stroke(ORANGE, width=4)
        )
        
        # 9. 두 막대 사이에 비교 표시
        #comparison_text = Text("Same Position!", font_size=14, color=WHITE)
        #comparison_text.next_to(target_bar, UP, buff=0.5)
        
        #self.play(Write(comparison_text))
        
        # 10. 두 막대를 번갈아 강조하여 비교 효과
        for _ in range(3):
            self.play(
                fifth_bar.animate.set_fill(opacity=1),
                target_bar.animate.set_fill(opacity=0.5),
                run_time=0.5
            )
            self.play(
                fifth_bar.animate.set_fill(opacity=0.5),
                target_bar.animate.set_fill(opacity=1),
                run_time=0.5
            )
        
        # 11. 최종 상태로 복원
        self.play(
            fifth_bar.animate.set_fill(YELLOW, opacity=0.8),
            target_bar.animate.set_fill(ORANGE, opacity=0.8)
        )
        
        self.wait(3)
