from manim import *
import numpy as np

config.frame_width = 12
config.pixel_width = 1280
config.pixel_height = 1024

class InterpolatingBars(Scene):
    def construct(self):
        # 10개 막대 초기 상태
        initial_10_values = [1/np.sqrt(i) for i in range(1, 11)]
        weight_10 = np.sum([1/i for i in range(1, 11)])
        initial_10_values = [i/np.sqrt(weight_10) for i in initial_10_values]
        
        # 20개 막대 최종 상태
        final_20_values = [1/np.sqrt(i) for i in range(1, 21)]
        weight_20 = np.sum([1/i for i in range(1, 21)])
        final_20_values = [i/np.sqrt(weight_20) for i in final_20_values]
        
        max_height = 2
        bar_width = 0.15
        initial_spacing = 0.3
        final_spacing = 0.2
        
        # 10개 막대 생성 (빨간색으로 시작)
        bars_10 = VGroup()
        labels_10 = VGroup()
        
        for i in range(10):
            height = initial_10_values[i] * max_height / 0.5
            bar = Rectangle(
                width=bar_width,
                height=height,
                fill_color=RED,
                fill_opacity=0.7,
                stroke_width=2
            )
            bar.shift(RIGHT * (i - 4.5) * initial_spacing)
            bar.shift(UP * height/2)
            
            label = Text(f"j_{i+1}", font_size=12)
            label.next_to(bar, DOWN, buff=0.1)
            
            bars_10.add(bar)
            labels_10.add(label)
        
        title1 = Text("Embezzling State (n=10)", font_size=18)
        title1.to_edge(UP)
        
        title2 = Text("Expanded Embezzling State (n=20)", font_size=18)
        title2.to_edge(UP)
        
        # 초기 상태 표시
        self.play(Write(title1))
        self.play(Create(bars_10))
        self.play(Write(labels_10))
        self.wait(2)
        
        # 제목 변경
        self.play(Transform(title1, title2))
        
        # 1단계: 기존 막대들의 간격을 넓혀서 공간 만들기
        expanded_positions = []
        for i in range(10):
            new_x = (i * 2 - 9) * final_spacing  # 홀수 인덱스 위치로 이동
            expanded_positions.append(RIGHT * new_x)
        
        # 기존 막대들을 새로운 위치로 이동 (간격 넓히기)
        move_animations = []
        for i, bar in enumerate(bars_10):
            target_pos = expanded_positions[i] + UP * (initial_10_values[i] * max_height) / 2
            move_animations.append(bar.animate.move_to(target_pos))
        
        self.play(AnimationGroup(*move_animations), run_time=1.5)
        
        # 2단계: 사이사이에 새로운 막대들 생성
        new_bars = VGroup()
        new_labels = VGroup()
        
        # 기존 막대 사이 + 양 끝에 새로운 막대 배치
        new_indices = []
        
        # 첫 번째 막대 앞에 하나
        new_indices.append(0)
        # 기존 막대들 사이에 하나씩
        for i in range(9):
            new_indices.append((i + 1) * 2)
        
        for idx_in_new, final_idx in enumerate(new_indices):
            if final_idx < 20:  # 20개 범위 내에서만
                height = final_20_values[final_idx] * max_height / 0.5
                bar = Rectangle(
                    width=bar_width,
                    height=height,
                    fill_color=RED,
                    fill_opacity=0.7,
                    stroke_width=2
                )
                # 사이사이 위치 (짝수 인덱스)
                x_pos = (final_idx - 9.5) * final_spacing
                bar.move_to(RIGHT * x_pos + UP * height/2)
                
                label = Text(f"j_{final_idx+1}", font_size=10)
                label.next_to(bar, DOWN, buff=0.1)
                
                new_bars.add(bar)
                new_labels.add(label)
        
        # 3단계: 새로운 막대들을 사이사이에서 나타나게 하기
        self.play(
            AnimationGroup(
                *[FadeIn(bar, shift=UP*0.5) for bar in new_bars],
                *[FadeIn(label) for label in new_labels]
            ),
            run_time=2
        )
        
        # 4단계: 기존 라벨들 업데이트
        updated_labels = VGroup()
        label_updates = []
        
        for i, label in enumerate(labels_10):
            new_text = f"j_{i*2+1}"  # 홀수 번째로 변경
            new_label = Text(new_text, font_size=10)
            new_label.move_to(label.get_center())
            label_updates.append(Transform(label, new_label))
        
        self.play(AnimationGroup(*label_updates))
        
        # 5단계: 모든 막대를 순차적으로 강조하여 최종 분포 보여주기
        all_bars = VGroup(*bars_10, *new_bars)
        for bar in all_bars:
            self.play(
                bar.animate.set_fill(opacity=1),
                rate_func=there_and_back,
                run_time=0.15
            )
        
        self.wait(3)
