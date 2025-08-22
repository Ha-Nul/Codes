from manim import *
import numpy as np

config.frame_width = 12
config.pixel_width = 1280
config.pixel_height = 1024

class SquaredMajorizationVisualization(Scene):
    def construct(self):
        # 초기 값들 (maximally entangled)
        initial_values = [1/np.sqrt(10) for _ in range(10)] 
        # 변화된 값들 (embezzling state)
        changed_values = [1/np.sqrt(i) for i in range(1,11)]
        weight = np.sum([1/i for i in range(1,11)])
        changed_values = [i/np.sqrt(weight) for i in changed_values]
        
        # **수정된 부분**: 각 원소의 제곱 계산
        initial_values_squared = np.array(initial_values)**2
        changed_values_squared = np.array(changed_values)**2
        
        # **수정된 부분**: 제곱의 누적합 계산
        max_entangled_cumsum_squared = np.cumsum(initial_values_squared)
        embezzling_cumsum_squared = np.cumsum(changed_values_squared)
        
        # 전체 최대값 구하기 (개별 제곱값과 제곱의 누적합 포함)
        all_values = (list(initial_values_squared) + list(changed_values_squared) + 
                     list(max_entangled_cumsum_squared) + list(embezzling_cumsum_squared))
        global_max = max(all_values)
        
        # 동일한 스케일링 적용
        uniform_scale = 3.0 / global_max
        
        max_height = 1
        bar_width = 0.15
        spacing = 0.25
        
        bars = VGroup()
        labels = VGroup()
        
        # **수정된 부분**: 초기 막대들을 제곱값으로 생성
        for i in range(10):
            height = initial_values_squared[i] * uniform_scale * max_height
            bar = Rectangle(
                width=bar_width,
                height=height,
                fill_color=BLUE,
                fill_opacity=0.7,
                stroke_width=2
            )
            bar.shift(RIGHT * (i - 4.5) * spacing)
            bar.shift(UP * height/2)
            
            #label = Text(f"j_{i+1}", font_size=12)
            #label.next_to(bar, DOWN, buff=0.1)
            
            bars.add(bar)
            #labels.add(label)
        
        #title1 = Text("Maximally Entangled State (Squared Values)", font_size=16)
        #title1.to_edge(UP)
        
        #self.play(Write(title1))
        self.play(Create(bars))
        self.play(Write(labels))
        self.wait(2)
        
        # **수정된 부분**: Embezzling State의 제곱값으로 변환
        #title2 = Text("Embezzling State (Squared Values)", font_size=16)
        #title2.to_edge(UP)
        #self.play(Transform(title1, title2))
        
        for i, bar in enumerate(bars):
            new_height = changed_values_squared[i] * uniform_scale * max_height
            new_bar = Rectangle(
                width=bar_width,
                height=new_height,
                fill_color=RED,
                fill_opacity=0.7,
                stroke_width=2
            )
            new_bar.shift(RIGHT * (i - 4.5) * spacing)
            new_bar.shift(UP * new_height/2)
            
            self.play(Transform(bar, new_bar), run_time=0.3)
        
        self.wait(2)
        
        # **수정된 부분**: 제곱의 누적합으로 Majorization 시각화
        #title3 = Text("Majorization: Cumulative Sum of Squared Values", font_size=16)
        #title3.to_edge(UP)
        #self.play(Transform(title1, title3))
        
        # Embezzling 제곱의 누적합으로 변환
        for i, bar in enumerate(bars):
            cumsum_squared_value = embezzling_cumsum_squared[i]
            new_height = cumsum_squared_value * uniform_scale * max_height
            
            new_bar = Rectangle(
                width=bar_width,
                height=new_height,
                fill_color=PURPLE,
                fill_opacity=0.7,
                stroke_width=2
            )
            new_bar.shift(RIGHT * (i - 4.5) * spacing)
            new_bar.shift(UP * new_height/2)
            
            self.play(Transform(bar, new_bar), run_time=0.4)
        
        self.wait(2)
        
        # **수정된 부분**: Maximally entangled 제곱의 누적합
        comparison_bars = VGroup()
        
        for i in range(10):
            cumsum_squared_value = max_entangled_cumsum_squared[i]
            height = cumsum_squared_value * uniform_scale * max_height
            
            bar_comp = Rectangle(
                width=bar_width,
                height=height,
                fill_color=BLUE,
                fill_opacity=0.5,
                stroke_width=2
            )
            bar_comp.shift(RIGHT * (i - 4.5) * spacing)
            bar_comp.shift(UP * height/2)
            
            comparison_bars.add(bar_comp)
        
        self.play(FadeIn(comparison_bars))
        
        # 확률 정보 표시
        prob_info = VGroup(
            Text("Squared Values = Probabilities", font_size=12),
            Text(f"Max Ent. Final Prob: {max_entangled_cumsum_squared[-1]:.3f}", font_size=11),
            Text(f"Embezzling Final Prob: {embezzling_cumsum_squared[-1]:.3f}", font_size=11),
            Text(f"Scale Factor: {uniform_scale:.3f}", font_size=11)
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        prob_info.to_corner(DL, buff=0.5)
        
        self.play(Write(prob_info))
        
        # 범례
        legend_blue = Rectangle(width=0.3, height=0.1, fill_color=BLUE, fill_opacity=0.5)
        legend_blue.to_corner(UR, buff=1).shift(LEFT * 2)
        legend_blue_text = Text("Max Entangled (Prob. Cumsum)", font_size=12)
        legend_blue_text.next_to(legend_blue, RIGHT, buff=0.1)
        
        legend_purple = Rectangle(width=0.3, height=0.1, fill_color=PURPLE, fill_opacity=0.7)
        legend_purple.next_to(legend_blue, DOWN, buff=0.2)
        legend_purple_text = Text("Embezzling (Prob. Cumsum)", font_size=12)
        legend_purple_text.next_to(legend_purple, RIGHT, buff=0.1)
        
        self.play(
            Create(legend_blue), Write(legend_blue_text),
            Create(legend_purple), Write(legend_purple_text)
        )
        
        self.wait(3)
