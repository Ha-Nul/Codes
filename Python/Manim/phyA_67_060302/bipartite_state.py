from manim import *
import numpy as np

# 코드 최상단에 추가
config.frame_width = 12  # 기본값 14.222
#config.frame_height = 6.75  # 16:9 비율 유지[6][14]
config.pixel_width = 1280  # 실제 픽셀 너비
config.pixel_height = 1024  # 실제 픽셀 높이[17]


class SimpleBipartiteState(Scene):
    def construct(self):
        # 수동으로 Rectangle들을 만들어서 막대 그래프 구현
        values = [0.316] * 10
        max_height = 2
        bar_width = 0.2
        spacing = 0.3
        
        bars = VGroup()
        labels = VGroup()
        
        for i in range(10):
            # 각 막대를 Rectangle로 생성
            height = values[i] * max_height / 0.4  # 0.4가 최대값
            bar = Rectangle(
                width=bar_width,
                height=height,
                fill_color=BLUE,
                fill_opacity=0.7,
                stroke_width=1
            )
            bar.shift(RIGHT * (i - 4.5) * spacing)
            #bar.shift(DOWN)  # 바닥에서 시작하도록
            
            # 숫자 라벨
            label = Text(str(i), font_size=14)
            label.next_to(bar, DOWN, buff=0.2)
            
            # 값 라벨
            value_label = Text("0.316", font_size=16)
            value_label.next_to(bar, UP, buff=0.1)
            
            bars.add(bar)
            labels.add(label)
            labels.add(value_label)

        
        # 제목
        title = Text("Bipartite Maximally Entangled State", font_size=24)
        title.to_edge(UP)
        
        # 애니메이션
        self.play(Write(title))
        self.play(Create(bars))
        self.play(Write(labels))
        self.wait(2)
