from manim import *
import numpy as np

config.frame_width = 12  # 기본값 14.222
#config.frame_height = 6.75  # 16:9 비율 유지[6][14]
config.pixel_width = 1280  # 실제 픽셀 너비
config.pixel_height = 1024  # 실제 픽셀 높이[17]

class ChangingBipartiteState(Scene):
    def construct(self):
        # 초기 값들 (maximally entangled)
        initial_values = [1/np.sqrt(10) for i in range(1,11,1)] 
        # 변화된 값들 (non-maximally entangled)
        changed_values = [1/np.sqrt(i) for i in range(1,11,1)]
        weight = np.sum([1/i for i in range(1,11,1)])
        changed_values = [i/np.sqrt(weight) for i in changed_values]
        
        max_height = 2
        bar_width = 0.2
        spacing = 0.3
        
        bars = VGroup()
        labels = VGroup()
        
        # 초기 막대들 생성
        for i in range(10):
            height = initial_values[i] * max_height / 0.5
            bar = Rectangle(
                width=bar_width,
                height=height,
                fill_color=BLUE,
                fill_opacity=0.7,
                stroke_width=2
            )
            bar.shift(RIGHT * (i - 4.5) * spacing)
            bar.shift(UP * height/2)
            
            label = Text(r"j_" + str(i), font_size=12)
            label.next_to(bar, DOWN, buff=0.1)
            
            bars.add(bar)
            labels.add(label)
        
        title1 = Text("Maximally Entangled State", font_size=18)
        title1.to_edge(UP)
        
        # 첫 번째 상태 보여주기
        self.play(Write(title1))
        self.play(Create(bars))
        self.play(Write(labels))
        self.wait(2)
        
        # 제목 변경
        title2 = Text("Embezzling State", font_size=18)
        title2.to_edge(UP)
        self.play(Transform(title1, title2))
        
        # 막대 높이 변경 애니메이션
        for i, bar in enumerate(bars):
            new_height = changed_values[i] * max_height / 0.5
            new_bar = Rectangle(
                width=bar_width,
                height=new_height,
                fill_color=RED,
                fill_opacity=0.7,
                stroke_width=2
            )
            new_bar.shift(RIGHT * (i - 4.5) * spacing)
            new_bar.shift(UP * new_height/2)
            
            self.play(Transform(bar, new_bar), run_time=0.5)
        
        self.wait(3)
